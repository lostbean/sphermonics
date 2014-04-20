{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Texture.SH.RotationSO3 where

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

import           Texture.HyperSphere
import           Texture.Orientation

import           Texture.SH.Harmonics
import           Texture.SH.Pyramid
import           Texture.SH.SupportFunctions

import           Debug.Trace
dbg s x = trace (show s L.++ show x) x

-- ============================== Spherical Harmonics Rotation ===========================

-- | Active rotation matrix series for Spherical Harmonics.
vecActiveRotMatrixSH :: Int -> SO3 -> Vector (M.Matrix (Complex Double))
vecActiveRotMatrixSH = calcVecRotMatrixSH id

-- | Passive rotation matrix series for Spherical Harmonics.
vecPassiveRotMatrixSH :: Int -> SO3 -> Vector (M.Matrix (Complex Double))
vecPassiveRotMatrixSH = calcVecRotMatrixSH conjugate

-- | Active rotation matrix series for Real Spherical Harmonics.
vecActiveRotMatrixRealSH :: Int -> SO3 -> Vector (M.Matrix Double)
vecActiveRotMatrixRealSH = calcVecRotMatrixRealSH id

-- | Passive rotation matrix series for Real Spherical Harmonics.
vecPassiveRotMatrixRealSH :: Int -> SO3 -> Vector (M.Matrix Double)
vecPassiveRotMatrixRealSH = calcVecRotMatrixRealSH conjugate

-- | Rotation matrix series for spherical harmonics from 2l = 0 until 2l = 2n. More
-- info on PhD Thesis J.K. Mason eq (52).
calcVecRotMatrixSH :: (Complex Double -> Complex Double) -> Int -> SO3 ->
                      Vector (M.Matrix (Complex Double))
calcVecRotMatrixSH foo llen r = let
  z = genSHFunc (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = foo $ z %! (N (2*l), L lb, MF mu)
    in ((k1 * c1) :+ 0) * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps

  buildM l = let
    ms = [l, l-1 .. (-l)]
    ps = [(i, j) | i <- ms, j <- ms]
    in ((2*l+1) >< (2*l+1)) $ map (calcAtM l) ps
  in V.generate (llen+1) buildM

-- | Same as for 'calcVecRotMatrixSH' but for real Spherical Harmonics.
calcVecRotMatrixRealSH :: (Complex Double -> Complex Double) -> Int -> SO3 ->
                          Vector (M.Matrix Double)
calcVecRotMatrixRealSH foo l r = let
  vu = calcVecRotMatrixSH foo l r
  func i u = let
    c  = realPartMatrixSH i
    c' = M.trans $ M.mapMatrix conjugate c
    in M.mapMatrix realPart $ c' <> u <> c
  in V.imap func vu

-- ======================== Matrix multiplication on L Stack Level =======================

-- | Applies matrix multiplication for each MF column along the L stack of the 'Pyramid'.
multLStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (L, MF) a -> Pyramid (L, MF) a
multLStack getMs pyC = let
  l       = getMaxStack pyC
  ls      = [0 .. l]
  toVec   = U.fromList . M.toList
  func li = (getMs li) <> (getLSlice pyC (L li))
  pyra    = U.concat $ map (toVec . func) ls
  in pyC { pyramid = pyra }

getLSlice :: (U.Unbox a, Storable a)=> Pyramid (L, MF) a -> L -> M.Vector a
getLSlice pyC l@(L li) = toVec $ slicer $ pyramid pyC
  where
    toVec  = M.fromList . U.toList
    slicer = U.slice (getKeyPos (l, MF li)) (2 * li + 1)

-- ================================ Useful SH functions ==================================

-- | Applies active rotation on a Spherical Harmonics.
rotActiveSH :: SO3 -> Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) (Complex Double)
rotActiveSH rot = \sh -> let l = getMaxStack sh in multLStack (vecActiveRotMatrixSH l rot V.!) sh

-- | Applies passive rotation on a Spherical Harmonics.
rotPassiveSH :: SO3 -> Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) (Complex Double)
rotPassiveSH rot = \sh -> let l = getMaxStack sh in multLStack (vecPassiveRotMatrixSH l rot V.!) sh

-- | Applies active rotation on a Real Spherical Harmonics.
rotActiveRealSH :: SO3 -> Pyramid (L, MF) Double -> Pyramid (L, MF) Double
rotActiveRealSH rot = \sh -> let l = getMaxStack sh in multLStack (vecActiveRotMatrixRealSH l rot V.!) sh

-- | Applies passive rotation on a Real Spherical Harmonics.
rotPassiveRealSH :: SO3 -> Pyramid (L, MF) Double -> Pyramid (L, MF) Double
rotPassiveRealSH rot = \sh -> let l = getMaxStack sh in multLStack (vecPassiveRotMatrixRealSH l rot V.!) sh

-- ================================ Conversion to real Pyramid ===========================

-- | Transforms a complex Pyramid (L, MF) to its real equivalent.
fromComplexSH :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexSH pyC = generatePyramid (k . realPart . func) lmax
  where
    lmax = getMaxStack pyC
    k = (/ (sqrt 2))
    i = powerComplex
    func (l, mf)
      | mf > 0 = i (unL l)     * (s*pyC %! (l, mf) + (pyC %! (l, -mf)))
      | mf < 0 = i (unL $ l-1) * (s*pyC %! (l, mf) - (pyC %! (l, -mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1

-- =================== Conversion to real numbers transformation matrix ==================

realPartMatrixSH :: Int -> M.Matrix (Complex Double)
realPartMatrixSH l = let
  calcAtM (i, j)
    | i == 0 && j == 0 = 1 :+ 0
    | abs i /= abs j   = 0 :+ 0
    | j < 0 && i < 0   = 0 :+ (-s / (sqrt 2))
    | j < 0 && i > 0   = 0 :+ (1 / (sqrt 2))
    | j > 0 && i < 0   = (s / (sqrt 2)) :+ 0
    | otherwise        = (1 / (sqrt 2)) :+ 0  --  j > 0 && i > 0
    where s = if even i then 1 else (-1)
  in ((2*l+1) >< (2*l+1)) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]

-- ======================================== Test ========================================

-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
-- Real SH Rotation matrix are broken. fromComplexSH doesn't work for SH but fromComplexSH_wrong works
-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

{--
rotPassiveSH + (SO3 0 0 pi) => active rotation
rotActiveSH + (SO3 0 0 (\x -> 2pi - x)) => active rotation
--}

testRotSH :: IO ()
testRotSH = let
  g1  = SO2 (pi/4) (pi/2)
  g2  = SO2 (pi/3) (1.5*pi)
  rot = SO3 (pi/3) (pi/3) (pi)
  in do
    plotSHPoints [g1, g2] [rot] [rot]
    plotSH_C "initial" [g1, g2] fromComplexSH_wrong 
    plotSH   "initial" [g1, g2] id
    plotSH_C "active"  [g1, g2] (fromComplexSH_wrong . rotActiveSH rot)
    plotSH   "active"  [g1, g2] (rotActiveRealSH rot)
    plotSH_C "passive" [g1, g2] (fromComplexSH_wrong . rotPassiveSH rot)
    plotSH   "passive" [g1, g2] (rotPassiveRealSH rot)

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
  s = 2 * l + 1  -- same as (N+1)^2 for N=0,2,4...
  in (s >< s) $ map calcAtM [(i, j) | i <- [-l..l], j <- [-l..l]]

vecZRotMatrixHSH :: Int -> Double -> Vector (M.Matrix Double)
vecZRotMatrixHSH n omega = let
  lmax = n `quot` 2
  func l = rotMZReal l omega
  in V.generate lmax func

fromComplexSH_wrong :: Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double
fromComplexSH_wrong pyC = generatePyramid (realPart . func) lmax
  where
    lmax = getMaxStack pyC
    kr = (1/(sqrt 2)) :+ 0
    ki = 0 :+ (-1/(sqrt 2))
    func (l, mf)
      | mf > 0 = kr * (pyC %! (l, mf)  + s*(pyC %! (l, -mf)))
      | mf < 0 = ki * (pyC %! (l, -mf) - s*(pyC %! (l,  mf)))
      | otherwise = pyC %! (l, mf)
      where s = if even (unMF mf) then 1 else -1
