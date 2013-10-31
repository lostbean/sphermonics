{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Hammer.Texture.SH.Rotation where

{--
( rotM
       , rotMReal
       , rotMZ
       , rotMZReal
       , multLStack
       )
--}

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Packed         as M

import           Data.Complex
import           Foreign.Storable      (Storable)
import           Data.Vector           (Vector)
import           Data.Packed.Matrix    ((><))
import           Numeric.LinearAlgebra ((<>))
import           Numeric.Container     (Product)

-- import           Hammer.Math.Algebra
import           Hammer.Texture.SH.Pyramid
import           Hammer.Texture.SH.HyperSphere
import           Hammer.Texture.SH.SupportFunctions
import           Hammer.Texture.SphericalHarmonics

import           Debug.Trace
dbg s x = trace (show s L.++ show x) x


cM :: Int -> M.Matrix (Complex Double)
cM l = let
  calcAtM (i, j)
    | i == 0 && j == 0 = 1 :+ 0
    | abs i /= abs j   = 0 :+ 0
    | j < 0 && i < 0   = 0 :+ (-s / (sqrt 2))
    | j < 0 && i > 0   = 0 :+ (1 / (sqrt 2))
    | j > 0 && i < 0   = (s / (sqrt 2)) :+ 0
    | j > 0 && i > 0   = (1 / (sqrt 2)) :+ 0
    where s = if even i then 1 else (-1)
  in ((2*l+1) >< (2*l+1)) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]

cMVec :: Int -> Vector (M.Matrix (Complex Double))
cMVec l = V.generate (l+1) cM

rotMReal :: Int -> SO3 -> Vector (M.Matrix Double)
rotMReal l r = let
  vu = rotM l r
  func i u = let
    c  = cM i
    c' = M.trans $ M.mapMatrix conjugate c
    in M.mapMatrix realPart $ c' <> u <> c
  in V.imap func vu

-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotM :: Int -> SO3 -> Vector (M.Matrix (Complex Double))
rotM llen r = let
  z = genSHFunc (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = z %! (N (2*l), L lb, MF mu)
    in ((k1 * c1) :+ 0) * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps
  mL l = ((2*l+1) >< (2*l+1)) $ map (calcAtM l) [(i, j) | i <- [-l .. l], j <- [-l .. l]]
  in V.generate (llen+1) mL


multLStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (L, MF) a -> Pyramid (L, MF) a
multLStack getMs pyC = let
  l       = maxStack pyC
  ls      = [0 .. l]
  toVec   = U.fromList . M.toList
  func li = (getMs li) <> (getLSlice pyC (L li))
  pyra    = U.concat $ map (toVec . func) ls
  in pyC { pyramid = pyra }

multNStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
multNStack getMs pyC = let
  n       = maxStack pyC
  ns      = [0, 2 .. n]
  func ni = multLStack getMs (getNSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

getLSlice :: (U.Unbox a, Storable a)=> Pyramid (L, MF) a -> L -> M.Vector a
getLSlice pyC l@(L li) = toVec $ slicer $ pyramid pyC
  where
    toVec  = M.fromList . U.toList
    slicer = U.slice (getKeyPos (l, MF (-li))) (2 * li + 1)

getNSlice :: (U.Unbox a, Storable a)=> Pyramid (N, L, MF) a -> N -> Pyramid (L, MF) a
getNSlice pyC nt = generatePyramid func ni
  where
    vec = slicer $ pyramid pyC
    n@(N ni) = if even (unN nt) then nt else nt-1
    func lmf = vec U.! getKeyPos (lmf :: (L, MF))
    slicer = U.slice (getKeyPos (n, L 0, MF 0)) ((ni + 1)^2)

fromComplexNLMF :: Pyramid (N, L, MF) (Complex Double) -> Pyramid (N, L, MF) Double
fromComplexNLMF pyC = let
  n       = maxStack pyC
  ns      = [0, 2 .. n]
  func ni = fromComplexLMF2 (getNSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

fromComplexLMF :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexLMF pyC = generatePyramid (realPart . func) lmax
  where
    lmax = maxStack pyC
    kr = (1/(sqrt 2)) :+ 0
    ki = 0 :+ (-1/(sqrt 2))
    func (l, mf)
      | mf > 0 = kr * (pyC %! (l, mf)  + s*(pyC %! (l, -mf)))
      | mf < 0 = ki * (pyC %! (l, -mf) - s*(pyC %! (l,  mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1

fromComplexLMF2 :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexLMF2 pyC = generatePyramid (k . realPart . func) lmax
  where
    lmax = maxStack pyC
    k = (/ (sqrt 2))
    i = powerComplex
    func (l, mf)
      | mf > 0 = i (unL l)     * (s*pyC %! (l, mf) + (pyC %! (l, -mf)))
      | mf < 0 = i (unL $ l-1) * (s*pyC %! (l, mf) - (pyC %! (l, -mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1

-- ======================================== Test ========================================

testRotHSH :: IO ()
testRotHSH = let
  foo  = (rotM 18 (SO3 (pi/4) (pi/4) 0) V.!)
  func = fromComplexNLMF . multNStack foo
  in testHSH_C func

testHSH_C :: (Pyramid (N,L,MF) (Complex Double) -> Pyramid (N,L,MF) Double) -> IO ()
testHSH_C func = let
  xs :: Vector (Complex Double, SO3)
  xs  = V.fromList [(10 :+ 0, SO3 (pi/2) (pi/2) (pi/2)), (10 :+ 0, SO3 pi (pi/4) 0)]
  c   = findSHCoefWeight 18 xs
  vtk = renderSO3SolidVTK (evalSH $ func c)
  in writeQuater "Rot-SHS-test" vtk

testSHRotC :: Pyramid (L, MF) (Complex Double)
testSHRotC = let
  xs :: Vector (Complex Double, SO2)
  xs = V.fromList [(10 :+ 0, (SO2 0 0)), (10 :+ 0, SO2 (pi/4) (pi/4)), (10 :+ 0, SO2 (pi/2) (pi/2))]
  in findSHCoefWeight 10 xs

testSHRot :: Pyramid (L, MF) (Double)
testSHRot = let
  xs :: Vector (Double, SO2)
  xs = V.fromList [(10, SO2 0 0), (10, SO2 (pi/4) (pi/4)), (10, SO2 (pi/2) (pi/2))]
  in findSHCoefWeight 10 xs

-- ======================================== Trash ========================================

rotMZ :: Int -> Double -> M.Matrix (Complex Double)
rotMZ l omega = let
  d = M.fromList $ map (\m -> cis (fromIntegral (m) * omega)) [-l .. l]
  in M.diagRect (0 :+ 0) d (2*l+1) (2*l+1)

rotMZReal :: Int -> Double -> M.Matrix Double
rotMZReal l omega = let
  calcAtM (i,j)
    |   i  == 0 && j == 0 = 1
    |   i  == j           = cos w
    | (-i) == j && i < 0  = sin w
    | (-i) == j && i > 0  = -(sin w)
    | otherwise           = 0
    where w = (fromIntegral i)*omega
  in ((2*l+1) >< (2*l+1)) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]
