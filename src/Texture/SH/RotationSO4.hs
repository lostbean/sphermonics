{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Texture.SH.RotationSO4 where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U
import qualified Data.HashMap.Strict as HM
import qualified Data.Packed         as M

import           Data.Complex
import           Foreign.Storable      (Storable)
import           Data.Vector           (Vector)
import           Data.HashMap.Strict   (HashMap)
import           Data.Packed.Matrix    ((><), (@@>), cols, buildMatrix)
import           Numeric.LinearAlgebra ((<>))
import           Numeric.Container     (Product, add, ident, ctrans, multiply, cmap)

import           Texture.HyperSphere
import           Texture.Orientation

import           Texture.SH.Harmonics
import           Texture.SH.Pyramid
import           Texture.SH.SupportFunctions
import           Texture.SH.RotationSO3

--import           Debug.Trace
--dbg s x = trace (show s L.++ show x) x

-- =========================== Hyper Spherical Harmonics Rotation ========================

-- | Rotation matrix series for hyper spherical harmonics from 2l = 0 until 2l = 2n. More
-- info on PhD Thesis J.K. Mason eq (52).
vecRotMatrixHSH :: Int -> SO3 -> SO3 -> Vector (M.Matrix (Complex Double))
vecRotMatrixHSH n r p = let
  lmax = n `quot` 2
  -- For some reason p and r had to swap positions, it mans that the active rotation is
  -- given by vecPassiveRotMatrixSH and the passive rotation is given by vecActiveRotMatrixSH  
  us   = vecActiveRotMatrixSH  lmax p
  us'  = vecPassiveRotMatrixSH lmax r
  in V.izipWith (\l u u' -> calcRotMatrixHSH u u' (2*l)) us us'

vecIDRotMatrixHSH :: Int -> Vector (M.Matrix (Complex Double))
vecIDRotMatrixHSH n = let
  lmax = n `quot` 2
  func l = let
    s  = (2*l+1) * (2*l+1)
    in ident s
  in V.generate lmax func

-- | Rotation matrix for hyper spherical harmonics at 2l. More info on PhD Thesis J.K.
-- Mason eq (52).
rotMatrixHSH :: Int -> SO3 -> SO3 -> M.Matrix (Complex Double)
rotMatrixHSH n r p = (vecRotMatrixHSH n r p) V.! (n `quot` 2)

-- | Rotation matrix series for Real hyper spherical harmonics from 2l = 0 until 2l = 2n.
vecRotMatrixRealHSH :: Int -> SO3 -> SO3 -> Vector (M.Matrix Double)
vecRotMatrixRealHSH n r p = let
  l  = n `quot` 2
  vu = vecRotMatrixHSH (2*l) r p
  vt = V.generate (2*l+1) realPartMatrixHSH
  func i u = let
    t  = M.diagBlock $ V.toList $ V.take (2*i+1) vt
    in M.mapMatrix realPart $ (ctrans t) <> u <> t
  in V.imap func vu

-- | Calculate the rotation matrix for Hyper Spherical Harmonics (SO4) based on two SO3
-- rotation matrix, one active and one passive. More info, see PhD Thesis J.K. Mason eq (52)
calcRotMatrixHSH :: M.Matrix (Complex Double) -> M.Matrix (Complex Double) -> Int -> M.Matrix (Complex Double)
calcRotMatrixHSH u u' n = let
  l  = n `quot` 2

  calcAtM (row, col) = let
    p' = ps V.! row
    p  = ps V.! col
    us = [l, l-1 .. (-l)]
    in sum [ func p' p (m', m) (m''', m'')
           | m <- us, m' <- us, m'' <- us, m''' <- us
           ]

  func (lb', mu') (lb, mu) (m', m) (m''', m'') = let
    -- calcClebshGordan j1 m1 j2 m2 j m
    --c1 = getClebsh (l, m'' ) (l, m') (lb', mu')
    --c2 = getClebsh (l, m''') (l, m ) (lb , mu )
    c1 = calcClebshGordan (l, m'' ) (l, m') (lb', mu')
    c2 = calcClebshGordan (l, m''') (l, m ) (lb , mu )
    z  = u  @@> (l - m'  , l - m  )
    z' = u' @@> (l - m''', l - m'')
    in ((c1 * c2) :+ 0) * z * z'

  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]
  ms = [ (row, col) | col <- [0 .. s-1], row <- [0 .. s-1] ]
  s  = (2*l+1) * (2*l+1)
  in (s >< s) $ map calcAtM ms

-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

isNonZeroClebsh :: (Int, Int) -> (Int, Int) -> (Int, Int) -> Bool
isNonZeroClebsh (j1, m1) (j2, m2) (j, m) =
  m  == m1 + m2       &&
  j  >= abs (j1 - j2) &&
  j  <= abs (j1 + j2) &&
  j  >= abs m         &&
  j1 >= abs m1        &&
  j2 >= abs m2

getClebsh :: (Int, Int) -> (Int, Int) -> (Int, Int) -> Double
getClebsh (l1, m1) (l2, m2) (lb, mu)
  | isNonZeroClebsh (l1, m1) (l2, m2) (lb, mu) = maybe slow id fast
  | otherwise                                  = 0
  where
    slow = calcClebshGordan (l1, m1) (l2, m2) (lb, mu)
    fast = HM.lookup (l1, m1, l2, m2, lb, mu) memoClebsh

memoClebsh :: HashMap (Int, Int, Int, Int, Int, Int) Double
memoClebsh = HM.fromList useful
  where
    l = 6
    useful =
      [ ((l1, m1, l2, m2, lb, mu)
      , calcClebshGordan (l1, m1) (l2, m2) (lb, mu))
      | l1 <- [0     .. l  ]
      , m1 <- [(-l)  .. l  ]
      , l2 <- [0     .. l  ]
      , m2 <- [(-l)  .. l  ]
      , lb <- [0     .. 2*l]
      , mu <- [(-lb) .. lb ]
      , isNonZeroClebsh (l1, m1) (l2, m2) (lb, mu)
      ]

-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-- ======================== Matrix multiplication on N Stack Level =======================

-- | Applies matrix multiplication for each L x MF column along the N stack of the 'Pyramid'.
multNStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
multNStack getMs pyC = let
  n       = getMaxStack pyC
  ls      = [0 .. (n `quot` 2)]
  toVec   = U.fromList . M.toList
  func li = (getMs li) <> (getNSlice pyC (N $ 2*li))
  pyra    = U.concat $ map (toVec . func) ls
  in pyC { pyramid = pyra }

getNSlice :: (U.Unbox a, Storable a)=> Pyramid (N, L, MF) a -> N -> M.Vector a
getNSlice pyC n = toVec $ slicer $ pyramid pyC
  where
    l      = (unN n) `quot` 2
    toVec  = M.fromList . U.toList
    slicer = U.slice (getKeyPos (N (2*l), L 0, MF 0)) ((2*l+1)*(2*l+1))

-- =================== Matrix multiplication on N over L Stack Level =====================

-- | Applies matrix multiplication for each MF column along the N and L stacks of the 'Pyramid'.
multNLStack :: (U.Unbox a, Storable a, Product a)=> (Int -> M.Matrix a) ->
             Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
multNLStack getMs pyC = let
  n       = getMaxStack pyC
  ns      = [0, 2 .. n]
  func ni = multLStack getMs (getNLSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

getNLSlice :: (U.Unbox a, Storable a)=> Pyramid (N, L, MF) a -> N -> Pyramid (L, MF) a
getNLSlice pyC nt = generatePyramid func ni
  where
    vec = slicer $ pyramid pyC
    n@(N ni) = if even (unN nt) then nt else nt-1
    func lmf = vec U.! getKeyPos (lmf :: (L, MF))
    slicer = U.slice (getKeyPos (n, L 0, MF 0)) ((ni + 1)^(2 :: Int))

-- =================== Conversion to real numbers transformation matrix ==================

realPartMatrixHSH :: Int -> M.Matrix (Complex Double)
realPartMatrixHSH l = let
  calcAtM (i, j)
    | mi == 0 && mj == 0 = (ip l)
    | abs mi /= abs mj   = 0 :+ 0
    | mj > 0 && mi > 0   = (ip l)     * k * s
    | mj > 0 && mi < 0   = (ip l)     * k
    | mj < 0 && mi > 0   = (ip $ l-1) * k * s
    | otherwise          = (ip $ l-1) * (-k)
    where
      -- mi and mj range [l .. -l]
      mi = -i + l
      mj = -j + l
      s  = if even mi then 1 else (-1)
  ip = negativeImaginaryPower
  k  = 1 / (sqrt 2)
  ms = 2 * l + 1
  in buildMatrix ms ms calcAtM

toRealMatrix :: M.Matrix (Complex Double) -> M.Matrix Double
toRealMatrix m = let
  nsqrt = round $ sqrt $ (fromIntegral $ cols m :: Double)
  l = nsqrt - 1
  u = M.diagBlock $ map realPartMatrixHSH [0 .. l]
  in cmap realPart $ (ctrans u) `multiply` m `multiply` u

-- ================================ Conversion to real Pyramid ===========================

-- | Transforms a complex Pyramid (L, MF) to its real equivalent.
fromComplexHSH :: Pyramid (N, L, MF) (Complex Double) -> Pyramid (N, L, MF) Double
fromComplexHSH pyC = generatePyramid (realPart . func) lmax
  where
    lmax = getMaxStack pyC
    kr   = (1 / sqrt 2) :+ 0
    i    = imaginaryPower . unL
    func (n, l, mf)
      | mf > 0    = i  l      * kr * (s * pyC %! (n, l, mf)  + pyC %! (n, l, -mf))
      | mf < 0    = i (l - 1) * kr * (s * pyC %! (n, l, -mf) - pyC %! (n, l,  mf))
      | otherwise = i l       * pyC %! (n, l, mf)
      where s = if even (unMF mf) then 1 else -1

-- ================================ Useful HSH functions =================================

-- | Applies both active and passive rotations on a Hyper Spherical Harmonics. The active
-- one is a rotation on the crystal reference frame and the passive is on the sample frame.
rotRHSH :: (SO3, SO3) -> Pyramid (N, L, MF) Double -> Pyramid (N, L, MF) Double
rotRHSH (ac, pa) = \hsh -> let
  n      = getMaxStack hsh
  vecmat = V.map toRealMatrix $ vecRotMatrixHSH n ac pa
  in applyRotMatrixHSH vecmat hsh

-- | Applies both active and passive rotations on a Hyper Spherical Harmonics. The active
-- one is a rotation on the crystal reference frame and the passive is on the sample frame.
rotHSH :: (SO3, SO3) -> Pyramid (N, L, MF) (Complex Double) -> Pyramid (N, L, MF) (Complex Double)
rotHSH (ac, pa) = \hsh -> let n = getMaxStack hsh in applyRotMatrixHSH (vecRotMatrixHSH n ac pa) hsh

-- | Calculates the matrix that enforces symmetry on Hyper Spherical Harmonics.
symmRotMatrixHSH :: [(SO3, SO3)] -> Int -> Vector (M.Matrix (Complex Double))
symmRotMatrixHSH []   l = vecIDRotMatrixHSH  l
symmRotMatrixHSH symm l = let
  rot0:rots = map (\(ac, pa) -> vecRotMatrixHSH l ac pa) symm
  in L.foldl' (V.zipWith add) rot0 rots

-- | Calculates the matrix that enforces symmetry on Hyper Spherical Harmonics.
symmRotMatrixRHSH :: [(SO3, SO3)] -> Int -> Vector (M.Matrix Double)
symmRotMatrixRHSH symm = V.map toRealMatrix . symmRotMatrixHSH symm

-- | Applies a given rotation matrix to Hyper Spherical Harmonics.
applyRotMatrixHSH :: (Product a, U.Unbox a)=> Vector (M.Matrix a) -> Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
applyRotMatrixHSH rotvec = multNStack (rotvec V.!)

-- ======================================== Test ========================================

testRotHSH :: IO ()
testRotHSH = let
  g1   = SO3 (pi/2) (pi/2) (pi/2)
  g2   = SO3 (pi/4) (pi/4)  0
  rot  = SO3 (-pi/2) (pi/2)  0
  in do
    plotHSHPoints [g1, g2] [rot] [rot]
    plotHSH_C "initial" [g1, g2] fromComplexHSH
    plotHSH   "initial" [g1, g2] id
    plotHSH_C "active"  [g1, g2] (fromComplexHSH . rotHSH (rot, SO3 0 0 0))
    plotHSH   "active"  [g1, g2] (rotRHSH (rot, SO3 0 0 0))
    plotHSH_C "passive" [g1, g2] (fromComplexHSH . rotHSH (SO3 0 0 0, rot))
    plotHSH   "passive" [g1, g2] (rotRHSH (SO3 0 0 0, rot))

testSymmHSH :: IO ()
testSymmHSH = let
  g  = SO3 (0.77*pi) (0.3*pi) (0.4*pi)
  r1 = (SO3 0 0 0, SO3 (-pi/2) (pi/2) 0)
  r2 = (SO3 (-pi/2) (pi/2) 0, SO3 0 0 0)
  r3 = (SO3 0 0 0, SO3 0 0 0)
  symm = [r1, r2, r3]
  rot  = symmRotMatrixHSH symm 6
  (as, ps) = unzip symm
  in do
    plotHSHPoints [g] as ps
    plotHSH_C "symmetry"  [g] (fromComplexHSH . applyRotMatrixHSH rot)

-- ======================================== Trash ========================================

testRowCol0 :: N -> SO3 -> SO3 -> Double
testRowCol0 n r p = let
  rp  = quaternionToSO3 ((so3ToQuaternion r) #<= (so3ToQuaternion p))
  pr  = quaternionToSO3 ((so3ToQuaternion p) #<= (so3ToQuaternion r))
  vC1 = rotHSHCol0 n pr
  vC2 = head $ M.toColumns $ rotMatrixHSH (unN n) r p
  vR1 = rotHSHRow0 n rp
  vR2 = head $ M.toRows    $ rotMatrixHSH (unN n) r p
  tc  = V.imap (\i x -> abs $ magnitude $ x - (vC2 M.@> i) ) vC1
  tr  = V.imap (\i x -> abs $ magnitude $ x - (vR2 M.@> i) ) vR1
  in max (V.maximum tc) (V.maximum tr)

rotHSHRow0 :: N -> SO3 -> Vector (Complex Double)
rotHSHRow0 n r = let
  l = (unN n) `quot` 2
  z = genSHFunc (2*l) r
  func (lb, mu) = let
    k1 = sqrt 2 * pi / (fromIntegral $ 2*l+1)
    z1 = conjugate $ z %! (N (2*l), L lb, MF mu)
    in (k1 :+ 0) * z1
  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]
  in V.map func ps

rotHSHCol0 :: N -> SO3 -> Vector (Complex Double)
rotHSHCol0 n r = let
  l = (unN n) `quot` 2
  z = genSHFunc (2*l) r
  func (lb, mu) = let
    k1 = (-1)^lb * sqrt 2 * pi / (fromIntegral $ 2*l+1)
    z1 = z %! (N (2*l), L lb, MF mu)
    in (k1 :+ 0) * z1
  ps = V.fromList [ (lb, mu) | lb <- [0 .. 2*l], mu <- [lb, lb-1 .. (-lb)] ]
  in V.map func ps
