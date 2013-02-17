{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeSynonymInstances #-} 
{-# LANGUAGE FlexibleInstances #-} 


module Hammer.Texture.Harmonics.BaseFunctions
       where

import Data.Complex
import Foreign.Storable (Storable)

import Data.Packed.Vector as HV
import Data.Packed.Matrix as HM

import qualified Numeric.Container  as NC
import qualified Data.Vector        as V

import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x

-- | This module calculate base functions and coefficients used
-- in in ordinary or generalized  spherical or hyperspherical harmonics.
-- The equations are based on the following papers:
-- -  "On the efficient calculation of ordinary and generalized spherical harmonics"

newtype M = M Int deriving (Num, Eq, Ord, Enum, Show)
newtype N = N Int deriving (Num, Eq, Ord, Enum, Show)
newtype L = L Int deriving (Num, Eq, Ord, Enum, Show)

--type CMatrix = Array U DIM2 (Complex Double)
--type RMatrix = Array U DIM2 Double
--type RVector = Array U DIM1 Double


-- ========================= Calculate Wigner d Matrix ==================
-- TODO Fix for phi=0,2*pi,..2*n*pi and phi=pi/2..(2*n*pi)+pi/2
wigner :: M -> L -> Double -> Vector Double
wigner mm ll@(L l) phi = HV.fromList $ go (N $ l-1) [] --R.fromListUnboxed (R.ix1 (l+1))
  where
    go nn [] = let
      ls = [wignerInit mm ll phi]
      in go nn ls

    go nn@(N n) ls@[l1] = let
      l0 = nextWigner l1 0 mm nn ll phi
      in if n < 0 then ls else go (nn-1) (l0:ls)
    
    go nn@(N n) ls@(l1:l2:_) = let
      l0 = nextWigner l1 l2 mm nn ll phi
      in if n < 0 then ls else go (nn-1) (l0:ls)

-- | Find the a new element from the Wigner matrix recursively
-- based on n+1 and n+2 elements. Defined only for phi [0,pi].
nextWigner :: Double -> Double -> M -> N -> L -> Double -> Double
nextWigner w1 w2 (M m) (N n) ll phi = 
  let
    d  = alpha (n+1) ll 
    a2 = alpha (n+2) ll
    a1 = let
      m' = fromIntegral m
      n' = fromIntegral n
      in 2 * ((n' + 1) * cos phi - m') / sin phi
  in (a1 * w1 - a2 * w2) / d
{-# INLINE nextWigner #-}

wignerInit :: M -> L -> Double -> Double
wignerInit (M m) (L l) phi = let
  f1 = fromIntegral $ factorial (2*l)
  f2 = fromIntegral $ factorial (l+m) * factorial (l-m)
  f  = sqrt (f1 / f2)
  c  = (cos halfPhi)^(l+m) * (sin halfPhi)^(l-m)
  halfPhi = phi / 2
  --signal  = (-1)^(l-m)
  in f * c 

alpha :: Int -> L -> Double
alpha n (L l) = let
  q = fromIntegral $ (l+n)*(l-n+1)
  in sqrt q

-- | Generate a sqr matrix with size = 2l+1 where (0,0) -> (-l,-l) and
-- (size, size) -> (l,l). Expect L >= 0.
wignerMatrixWith :: (Element a)=> L -> Double -> (M -> N -> Double -> a) -> Matrix a
wignerMatrixWith ll@(L l) phi func = let
  rs    = fromRows [wigner (M i) ll phi | i <- [-l..l]]
  size  = 2*l+1
  sop (x,y) = (x-l, y-l)
  pos (x,y) = (x+l, y)
  symm ij
    | n < 0 && m < 0 = f $ rs @@> pos (-n,-m)
    | n < 0          = f $ (-1)^(abs $ n-m) * rs @@> pos (n,m)
    | otherwise      = f $ rs @@> pos (m,n)
    where (m, n) = sop ij
          f      = func (M m) (N n)
  in buildMatrix size size symm
     

-- ======================= Spherical Harmonics =============================

complexLegendreMatrix :: L -> Double -> Matrix (Complex Double)
complexLegendreMatrix l phi = let
  func (M m) (N n) = iexp (m-n)
  in wignerMatrixWith l phi func

realLegendreMatrix :: L -> Double -> Matrix Double
realLegendreMatrix l phi = wignerMatrixWith l phi (\_ _ x -> x)

bunge :: L -> Double -> Double -> Double -> Matrix (Complex Double)
bunge l phi1 phi phi2 = let
  func (M m) (N n) x = let
    e1 = eima m phi2
    p  = iexp (m-n) x
    e2 = eima n phi1
    in e1 * p * e2
  in wignerMatrixWith l phi func

getCoef :: L -> Double -> Double -> Double -> Matrix (Complex Double)
getCoef ll@(L l) phi1 phi phi2 = let
  k      = 2 * (fromIntegral l) + 1
  func x = (k :+ 0) * conjugate x
  in mapMatrix func (bunge ll phi1 phi phi2)

getRealCoef :: L -> Double -> Double -> Double -> Matrix Double
getRealCoef l phi1 phi phi2 = toRealMatrix $ getCoef l phi1 phi phi2

getRealCoefList :: L -> Double -> Double -> Double -> [Matrix Double]
getRealCoefList ll phi1 phi phi2= map (\l -> getRealCoef l phi1 phi phi2) [(L 0) .. ll]

eval :: Double -> Double -> Double -> Matrix (Complex Double) -> Complex Double
eval phi1 phi phi2 c = let
  ll@(L l) = getL c
  ijs = [(i,j) | i<-[0..2*l], j<-[0..2*l]]
  t   = bunge ll phi1 phi phi2
  func acc ij = acc + (c @@> ij * t @@> ij)
  in foldl func (0:+0) ijs

evalReal :: Double -> Double -> Double -> Matrix Double -> Double
evalReal phi1 phi phi2 c = let
  ll@(L l) = getL c
  ijs = [(i,j) | i<-[0..2*l], j<-[0..2*l]]
  t   = toRealMatrix $ bunge ll phi1 phi phi2
  func acc ij = acc + (c @@> ij * t @@> ij)
  in foldl func 0 ijs

evalRealList :: Double -> Double -> Double -> [Matrix Double] -> Double
evalRealList phi1 phi phi2 cs = let
  func acc x = acc + evalReal phi1 phi phi2 x
  in foldl func 0 cs

-- ======================= Support functions =========================== 
cplxSqrt :: (Floating a, Ord a)=> a -> Complex a
cplxSqrt q
  | q >= 0    = (sqrt q) :+ 0
  | otherwise = 0 :+ sqrt (-q)

fact :: Integer -> Integer
fact n = product [1..n]

iexp :: (Num a)=> Int -> a -> Complex a
iexp n q
  | i == 0 =   q  :+   0
  | i == 1 =   0  :+   q
  | i == 2 = (-q) :+   0
  | i == 3 =   0  :+ (-q)
  | otherwise = error "Non-exhaustive on function 'iexp'"
  where i = n `mod` 4
{-# INLINE iexp #-}
        
eima :: (Floating a)=> Int -> a -> Complex a
eima m a = let
  x = a * fromIntegral m
  in cos x :+ sin x

getL :: Matrix a -> L
getL matrix
  | cs == rs && m == 0 && d >= 0 = L d 
  | otherwise = error "[Hammer] Improper Spherical Harmonics matrix."
  where
    (d,m) = (cs - 1) `divMod` 2
    cs    = cols matrix
    rs    = rows matrix
    
-- | According to Bunge (eq. 14.36 and 14.37)
toRealMatrix :: Matrix (Complex Double) -> Matrix Double
toRealMatrix matrix = let
  k         = 1 / sqrt 2
  size      = 2 * l + 1
  (L l)     = getL matrix
  sop (x,y) = (x-l, y-l)
  pos (x,y) = (x+l, y+l)
  --signal m  = (-1)^(m `mod` 2)
  getij     = (matrix @@>) . pos
  func ij
    | m == 0 && n == 0 = getij (m,n)
    | m > 0     = iexp (m-n)   k * (getij (m,n) + getij (-m,-n))
    | otherwise = iexp (1+m-n) k * (getij (m,n) - getij (-m,-n))
    where (m,n) = sop ij
  in buildMatrix size size (realPart . func)

-- =========================== Plot ODF ===========================================
saveToGNUPlot :: (Show a, Storable a)=> FilePath -> Matrix a -> IO ()
saveToGNUPlot name m = let
  rs  = rows m
  cs  = cols m
  ijs = [[(i,j) | i <- [0..rs-1]] | j <- [1..cs-1]]
  ls  = map (map (\ij@(i,j) -> show i ++ " " ++ show j ++ " " ++ show (m @@> ij))) ijs
  out = unlines $ map unlines ls
  in writeFile name out 

saveToVTKPlot :: (Show a, Storable a)=> FilePath -> [Matrix a] -> IO ()
saveToVTKPlot name ms = let
  ppMatrix m = let 
    rs  = rows m
    cs  = cols m
    ijs = [(i,j) | i <- [0..rs-1], j <- [0..cs-1]]
    in unlines $ map (show . (m @@>)) ijs
  out = unlines $ map ppMatrix ms
  head = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> <VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\"> <ImageData WholeExtent=\"0 17 0 17 0 17\" Origin=\"0 0 0\" Spacing=\"5 5 5\"> <Piece Extent=\"0 17 0 17 0 17\"> <PointData Scalars=\"cell_scalars\"> <DataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"cell_scalars\" format=\"ascii\">"
  tail = "</DataArray> </PointData> </Piece> </ImageData> </VTKFile>"
  in writeFile name (head ++ out ++ tail)

plotPhi2SectionReal :: [Matrix Double] -> Double -> Matrix Double
plotPhi2SectionReal c phi2 = let
  k = (pi/2)/18
  func (i,j) = evalRealList (k*fromIntegral i) (k*fromIntegral j) phi2 c
  in buildMatrix 18 18 func

plotPhi2Section :: Matrix (Complex Double) -> Double -> Matrix (Complex Double)
plotPhi2Section c phi2 = let
  k = (pi/2)/18
  func (i,j) = eval (k*fromIntegral i) (k*fromIntegral j) phi2 c
  in buildMatrix 18 18 func

-- ========================= Check Properties =====================================
instance Element Bool
checkMatrix :: L -> Double -> Matrix Bool
checkMatrix ll phi = let
  matrix = complexLegendreMatrix ll phi
  rs = rows matrix
  cs = cols matrix
  l = (cs-1) `div` 2
  sop (x,y) = (x-l, y-l)
  pos (x,y) = (x+l, y+l)
  func ij = let
    (m,n) = sop ij
    a     = matrix @@> ij
    c     = matrix @@> pos (-n,-m)
    ac = magnitude $ a - c
    in rs == cs && ac < 0.0001
  in buildMatrix rs cs func

prop_angle_composition :: L -> Double -> Matrix Bool
prop_angle_composition l phi = let
  a = complexLegendreMatrix l (phi*0.3)
  b = complexLegendreMatrix l (phi*0.7)
  c = complexLegendreMatrix l phi
  m = (NC.multiply a b) `NC.sub` c
  in mapMatrix ((< 1e-8) . magnitude) m
  
-- ================= Fast factorial calculation ===================================

factorial :: (Integral a)=> a -> Integer
factorial n
  | n >= 200  = factorial' n
  | otherwise = fastFList V.! (fromIntegral n)

fastFList :: V.Vector Integer
fastFList = V.generate 200 factorial'

-- | Source: "[Haskell-cafe] faster factorial function via FFI?"
primePowerOf :: Integer -> Integer -> Integer
primePowerOf n p = (n - func p n) `div` (p - 1)
  where
    func _ 0 = 0
    func x i = let (q,r) = i `divMod` x in func x q + r

factorisedFactorial :: Integer -> [(Integer,Integer)]
factorisedFactorial n = [ (p, primePowerOf n p) | p <- takeWhile (<= n) primes ]

factorial' :: (Integral a)=> a -> Integer
factorial' = product . map (uncurry (^)) . factorisedFactorial . fromIntegral

-- List of primes. This is reasonably fast but doesn't have to be,
-- since we only need it this once to generate the lookup table
primes :: [Integer]
primes = (2 : [x | x <- [3,5..], isPrime x])

-- Efficient test presupposes the existence of primes
-- This works because to determine whether p is prime you only need
-- to know the primes strictly less than p (except for 2 of course!)
isPrime :: Integer -> Bool
isPrime x = null divisors
  where
    divisors      = [y | y <- onlyUpToSqrtX primes, x `mod` y == 0]
    onlyUpToSqrtX = fst . span (<= sqrtX)
    sqrtX         = floor (sqrt (fromIntegral x))
-- ================================================================================