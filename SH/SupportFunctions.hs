{-# LANGUAGE BangPatterns #-}

module Hammer.Texture.SH.SupportFunctions
       ( -- * Associated Legendre polynomials (Recursion)
         genAssLegenPyramid
       , genAssLegenFullPyramid
         -- * Associated Legendre polynomials (Hypergeometric function)
       , genAssLegenPyramidSlow
       , calcAssLegenSlow
         -- * Gegenbauer coefficient
       , genGegenbauerPyramid
         -- * Clebsh-Gordan coefficient
       , calcClebshGordan
         -- * Factorials
       , fact
       , logFact
       ) where

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM

import           Control.Monad.ST    (ST)
import           Math.Gamma (gamma)

import           Hammer.Texture.SH.Pyramid

import Debug.Trace

-- ===================== Associated Legendre Polynomials by 2F1 ========================== 

genAssLegenPyramidSlow :: L -> Double -> Pyramid (L, MF) Double
genAssLegenPyramidSlow (L li) x = generatePyramid (\(l, m) -> calcAssLegenSlow l m x) li

-- | Works for L values up to 10.
calcAssLegenSlow :: L -> MF -> Double -> Double
calcAssLegenSlow pl@(L l) pm@(MF m) x
  | m > 0     = k * calcPWithNegativeM pl (M (-m)) x * exp (a - b)
  | otherwise = calcPWithNegativeM pl (M m) x
  where
    a = logFact (l+m)
    b = logFact (l-m)
    k = if even m then 1 else -1

calcPWithNegativeM :: L -> M -> Double -> Double
calcPWithNegativeM (L l) (M m) x
  | r == 0 && k < 0 = 0   -- Solve singularity when x = -1
  | otherwise       = r**k / gamma c * f
  where
    f = f21 a b c z
    r = (1+x) / (1-x)
    z = (1-x) / 2
    k = (fromIntegral m) / 2
    a = fromIntegral (-l)
    c = fromIntegral (1 - m)
    b = fromIntegral (l + 1)

f21 :: Double -> Double -> Double -> Double -> Double
f21 a b c x = go 0 1
  where
    go j cj
      | j > 10    = 1
      | otherwise = let
        cj1 = cj * ((a+j)*(b+j)/(c+j)) * (x/(j+1))
        in cj1 + go (j+1) cj1

-- ============================= Generate Associate Legendre =============================

genAssLegenFullPyramid :: L -> Double -> Pyramid (L, MF) Double
genAssLegenFullPyramid ll x = generatePyramid func (unL ll)
  where
    pyP = genAssLegenPyramid ll x
    func (l, mf)
      | mf >= 0   = p
      | otherwise = pMinvL lm p
      where
        lm = (l, mf2m mf)
        p  = pyP %! lm

genAssLegenPyramid :: L -> Double -> Pyramid (L, M) Double
genAssLegenPyramid l x = let
  kmax = getMaxKey (unL l) :: (L, M)
  -- reversed sequence l -> m (top -> bottom)
  ps = [ (li, mi)
       | li <- [0 .. l]
       , mi <- let m = M (unL li) in [m, m-1 .. 0]]
  foo :: ST s (UM.MVector s Double)
  foo = do
    v <- UM.new (getLinSize kmax)
    mapM_ (fillP v x) ps
    return v
  vec = U.create foo
  in unsafeMkPyramid (unL l) vec

fillP :: UM.MVector s Double -> Double -> (L, M) -> ST s ()
fillP v x lm@(pl, pm)
  | l == 0 && m == 0 = goP00
  | l == m           = goDiag
  | l == m+1         = goSubDiag
  | m == 0           = goCentralLine
  | otherwise        = goFillLayer   -- m > 0
  where
    l = unL pl
    m = unM pm
    pos    = getKeyPos lm
    goP00  = UM.write v pos 1
    goDiag = UM.write v pos (pLL pl x)
    goSubDiag = do
      let prevl = pl - 1
      pll <- UM.read v (getKeyPos (prevl, pm))
      UM.write v pos (pLLUp (prevl) x pll)
    goCentralLine =  do
      pa <- UM.read v (getKeyPos (pl-1, pm))
      pb <- UM.read v (getKeyPos (pl-2, pm))
      UM.write v pos (pLUpM lm x pa pb)
    goFillLayer = do
      pa <- UM.read v (getKeyPos (pl-1, pm+1))
      pb <- UM.read v (getKeyPos (pl-1, pm-1))
      UM.write v pos (pML lm x pa pb)

{-# INLINE pLL #-}
pLL :: L -> Double -> Double
pLL (L l) x = let
  a = 0.5 * (fromIntegral l) * log (1-x*x)
  b = logFact (2*l)
  c = (fromIntegral l) * (log 2) + logFact l
  foo = if even l then id else negate
  in foo $ exp (a + b - c)

{-# INLINE pLLUp #-}
pLLUp :: L -> Double -> Double -> Double
pLLUp (L l) x pll = (fromIntegral $ 2*l+1) * x * pll

{-# INLINE pLUpM #-}
pLUpM :: (L, M) -> Double -> Double -> Double -> Double
pLUpM (L l0, M m) x pl1 pl2 = (vala - valb) / k
  where
    l    = l0 - 1
    k    = fromIntegral $ l-m+1
    vala = (fromIntegral $ 2 * l + 1) * x * pl1
    valb = (fromIntegral $ l + m) * pl2

{-# INLINE pML #-}
pML :: (L, M) -> Double -> Double -> Double -> Double
pML (L l, M m) x pbu pbl = let
  a = -sqrt (1-x*x)
  b = (pbu + (fromIntegral $ (l+m-1) * (l+m)) * pbl)
  c = (fromIntegral $ 2*m)
  in a * b / c

{-# INLINE pMinvL #-}
pMinvL :: (L, M) -> Double -> Double
pMinvL (L l, M m) plm = let
  k = exp (logFact (l - m) - logFact (l + m))
  in if even m then k * plm else -k * plm

-- ========================= Clebsch-Gordan coefficients =================================

-- | Following the http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html
-- constraints.
calcClebshGordan :: (Int, Int) -> (Int, Int) -> (Int, Int) -> Double
calcClebshGordan (j1, m1) (j2, m2) (j, m)
  | m  == m1 + m2       &&
    j  >= abs (j1 - j2) &&
    j  <= abs (j1 + j2) &&
    j  >= abs m         &&
    j1 >= abs m1        &&
    j2 >= abs m2        = delta * (exp (0.5 * (p1 + p3 - p2))) * suma
  | otherwise = 0
  where
    mapfact = U.sum . U.map logFact
    delta = if m == m1 + m2 then 1 else 0
    p1 = log (fromIntegral (2 * j + 1)) +
         logFact (j1 + j2 - j)  +
         logFact (j  + j1 - j2) +
         logFact (j  + j2 - j1)
    p2 = logFact (j1 + j2 + j + 1)
    p3 = logFact (j  - m ) + logFact (j  + m ) +
         logFact (j1 - m1) + logFact (j1 + m1) +
         logFact (j2 - m2) + logFact (j2 + m2)
    fk1 k = U.fromList [ (j1 + j2 - j - k), (j1 - m1 - k), (j2 + m2 - k) ]
    fk2 k = U.fromList [ (j  - j2 + m1 + k) , (j  - j1 - m2 + k) ]

    maxK = U.minimum $ fk1 0
    minK = max 0 (negate $ U.minimum $ fk2 0)
    ks   = U.enumFromN minK (maxK - minK + 1)
    func k = let
      s = if even k then 1 else -1
      d = logFact k + mapfact (fk1 k) + mapfact (fk2 k)
      in s * (exp (-d))
    suma = U.sum $ U.map func ks

-- ============================== Gegenbauer polynomials =================================

genGegenbauerPyramid :: N -> Double -> Pyramid (N, L) Double
genGegenbauerPyramid n x = let
  kmax = getMaxKey (unN n) :: (N, L)
  ps = [(ni, li)
       | ni <- [0 .. n]
       , li <- let l = L (unN ni) in [l, l-1 .. 0]]
  foo :: ST s (UM.MVector s Double)
  foo = do
    v <- UM.new (getLinSize kmax)
    mapM_ (fillGG v x) ps
    return v
  vec = U.create foo
  in unsafeMkPyramid 0 vec

fillGG :: UM.MVector s Double -> Double -> (N, L) -> ST s ()
fillGG v x nl@(n, l)
  | i == 0    = goP0
  | i == 1    = goP1
  | otherwise = go
  where
    a   = fromIntegral $ unL l + 1
    i   = unN n - unL l
    pos = getKeyPos nl
    func c1 c2 = let
      g1 = 2 * x * (i' + a - 1) * c1
      g2 = (i' + 2 * a - 2) * c2
      i' = fromIntegral i
      in (g1 - g2) / i'
    goP0 = UM.write v pos 1
    goP1 = UM.write v pos (2 * a * x)
    go = do
      c1 <- UM.read v (getKeyPos (n-1, l))
      c2 <- UM.read v (getKeyPos (n-2, l))
      UM.write v pos (func c1 c2)

calcGegenbauer :: Int -> Double -> Double -> Double
calcGegenbauer n0 a x = goG n0
  where
    goG !n
      | n <= 0 = 1
      | n == 1 = 2 * a * x
      | otherwise = let
        g1  = 2 * x * (nd + a - 1) * goG (n - 1)
        g2  = (nd + 2 * a - 2) * goG (n - 2)
        nd = fromIntegral n
        in (g1 - g2) / nd

testGegen :: Int -> Double -> Pyramid (N, L) Double
testGegen nmax x = let
  pG = genGegenbauerPyramid (N nmax) x
  func nl@(N n, L l) = pG %! nl - calcGegenbauer (n - l) (fromIntegral $ l + 1) x
  in generatePyramid func nmax

-- ================================= Factorial tables ====================================

factTable :: U.Vector Int
factTable = U.scanl' (*) 1 (U.enumFromN 1 10000)

fact :: Int -> Int
fact = (factTable U.!)

logFactTable :: U.Vector Double
logFactTable = U.scanl' (\acc i -> acc + log i) 0 (U.enumFromN 1 10000)

logFact :: Int -> Double
logFact = (logFactTable U.!)
