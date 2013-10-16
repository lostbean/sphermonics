{-# LANGUAGE BangPatterns #-}

module Hammer.Texture.SH.SupportFunctions
       ( -- * Associated Legendre polynomials (Recursion)
         genAssLegenPyramid
         -- * Associated Legendre polynomials (Hypergeometric function)
       , genAssLegenPyramidSlow
       , calcAssLegenSlow
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

genAssLegenPyramidSlow :: L -> Double -> Pyramid Double
genAssLegenPyramidSlow lm x = let
  struct = mkPyramidStruct lm HalfRange NoRange
  in generatePyramid (\(l, m, _) -> calcAssLegenSlow l m x) struct

-- | Works for L values up to 10.
calcAssLegenSlow :: L -> M -> Double -> Double
calcAssLegenSlow pl@(L l) pm@(M m) x
  | m > 0     = k * calcPWithNegativeM pl (-pm) x * a / b
  | otherwise = calcPWithNegativeM pl pm x
  where
    a = fromIntegral $ fact (l+m)
    b = fromIntegral $ fact (l-m)
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

genAssLegenPyramid :: L -> Double -> Pyramid Double
genAssLegenPyramid l x = let
  s  = mkPyramidStruct l HalfRange NoRange
  -- reversed sequence l -> m (top -> bottom)
  ps = [(li, mi, 0)
       | li <- [0 .. lmax s]
       , mi <- let (mb, mu) = getMRange s li in [mu, mu-1 .. mb]]
  foo :: ST s (UM.MVector s Double)
  foo = do
    v <- UM.new (linSize s)
    mapM_ (fillP s v x) ps
    return v
  vec = U.create foo
  in unsafeMkPyramid s vec

fillP :: PyramidStructure -> UM.MVector s Double -> Double -> LMN -> ST s ()
fillP s v x lmn@(pl, pm, pn)
  | l == 0 && m == 0 = goP00
  | l == m           = goDiag
  | l == m+1         = goSubDiag
  | m == 0           = goCentralLine
  | m >  0           = goFillLayer
  | otherwise        = goInvert      -- m < 0
  where
    l = unL pl
    m = unM pm
    pos    = getPos s lmn
    goP00  = UM.write v pos 1
    goDiag = UM.write v pos (pLL pl x)
    goSubDiag = do
      let prevl = pl - 1
      pll <- UM.read v (getPos s (prevl, pm, pn))
      UM.write v pos (pLLUp (prevl) x pll)
    goCentralLine =  do
      pa <- UM.read v (getPos s (pl-1, pm, pn))
      pb <- UM.read v (getPos s (pl-2, pm, pn))
      UM.write v pos (pLUpM lmn x pa pb)
    goFillLayer = do
      pa <- UM.read v (getPos s (pl-1, pm+1, pn))
      pb <- UM.read v (getPos s (pl-1, pm-1, pn))
      UM.write v pos (pML lmn x pa pb)
    goInvert = do
      let invLMN = (pl, -pm, pn)
      p <- UM.read v (getPos s invLMN)
      UM.write v pos (pMinvL invLMN p)

pLL :: L -> Double -> Double
pLL (L l) x = let
  a = sqrt ((1-x*x)^l)
  b = fromIntegral $ fact (2*l)
  c = fromIntegral $ 2^l * fact l
  foo = if even l then id else negate
  in foo (a * b / c)

pLLUp :: L -> Double -> Double -> Double
pLLUp (L l) x pll = (fromIntegral $ 2*l+1) * x * pll

pLUpM :: LMN -> Double -> Double -> Double -> Double
pLUpM (L l0, M m, _) x pl1 pl2 = (vala - valb) / k
  where
    l    = l0 - 1
    k    = fromIntegral $ l-m+1
    vala = (fromIntegral $ 2*l + 1) * x * pl1
    valb = (fromIntegral $ l + m) * pl2

pML :: LMN -> Double -> Double -> Double -> Double
pML (L l, M m, _) x pbu pbl = let
  a = -sqrt (1-x*x)
  b = (pbu + (fromIntegral $ (l+m-1) * (l+m)) * pbl)
  c = (fromIntegral $ 2*m)
  in a * b / c

pMinvL :: LMN -> Double -> Double
pMinvL (L l, M m, _) plm = let
  k = exp (logFact (l - m) - logFact (l + m))
  in if even m then k * plm else -k * plm

-- ========================= Clebsch–Gordan coefficients =================================

calcClebshGordan :: Int -> Int -> Int -> Int -> Int -> Int -> Double
calcClebshGordan j1 j2 m1 m2 j m
  | m == m1 + m2       &&
    abs (j1 - j2) <= j &&
    j <= j1 + j2       = delta * sqrt (p1 / p2) * sqrt p3 * suma
  | otherwise = 0
  where
    delta = if m == m1 + m2 then 1 else 0
    p1 = fromIntegral $ (2 * j + 1) *
         fact (j1 + j2 - j) *
         fact (j  + j1 - j2) *
         fact (j  + j2 - j1)

    p2 = fromIntegral $ fact (j1 + j2 + j + 1)

    p3 = fromIntegral $
         fact (j - m)   * fact (j + m) *
         fact (j1 - m1) * fact (j1 + m1) *
         fact (j2 - m2) * fact (j2 + m2)

    func k = let
      s = if even k then 1 else -1
      d = fact k *
          fact (j1 + j2 - j - k) *
          fact (j1 - m1 - k) *
          fact (j2 + m2 - k) *
          fact (j  - j2 + m1 + k) *
          fact (j  - j1 - m2 + k)
      in s / (fromIntegral d)
    minNonNeg = let
      ma = max (j - j2 + m1) (j - j1 - m2)
      in if ma > 0 then 0 else -ma
    maxNonNeg = let
      ma = min (j1 - m1) (j2 + m2)
      in min ma (j1 + j2 - j)
    nNonNeg = maxNonNeg - minNonNeg + 1
    ks = U.enumFromN minNonNeg nNonNeg
    suma = U.sum $ U.map func ks

-- ============================== Gegenbauer polynomials =================================

genGegenbauerPyramid :: LMN -> Double -> Pyramid Double
genGegenbauerPyramid (L l, M m, N n) x = let
  s  = mkPyramidStruct l HalfRange NoRange
  -- reversed sequence l -> m (top -> bottom)
  ps = [(li, 0, ni)
       | li <- [0 .. lmax s]
       , ni <- let (nb, nu) = getNRange s li in [nu, nu-1 .. nb]]
  foo :: ST s (UM.MVector s Double)
  foo = do
    v <- UM.new (linSize s)
    mapM_ (fillP s v x) ps
    return v
  vec = U.create foo
  in unsafeMkPyramid s vec

fillP :: PyramidStructure -> UM.MVector s Double -> Double -> LMN -> ST s ()
fillP s v x lmn@(pl, pm, pn)
  | l == 0 && m == 0 = goP00
  | otherwise        = goInvert      -- m < 0
  where
    l = unL pl
    m = unM pm
    a = pl + 1
    x = unN pn - unL pl
    pos    = getPos s lmn
    goP00  = UM.write v pos 1
    goDiag = UM.write v pos (pLL pl x)
    goSubDiag = do
      let prevl = pl - 1
      pll <- UM.read v (getPos s (prevl, pm, pn))
      UM.write v pos (pLLUp (prevl) x pll)

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


-- ================================= Factorial tables ====================================

factTable :: U.Vector Int
factTable = U.scanl' (*) 1 (U.enumFromN 1 10000)

fact :: Int -> Int
fact = (factTable U.!)

logFactTable :: U.Vector Double
logFactTable = U.scanl' (\acc i -> acc + log i) 0 (U.enumFromN 1 10000)

logFact :: Int -> Double
logFact = (logFactTable U.!)
