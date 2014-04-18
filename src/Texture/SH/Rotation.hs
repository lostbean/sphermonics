{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Texture.SH.Rotation
       ( -- * Spherical Harmonics Rotation
         rotActiveSH
       , rotPassiveSH
       , rotActiveRealSH
       , rotPassiveRealSH
         -- * Hyper Spherical Harmonics Rotation
       , rotHSH
       , symmRotMatrixHSH
       , applyRotMatrixHSH
         -- * HSH Rotation Matrices
       , vecRotMatrixHSH
       , rotMatrixHSH
       , vecRotMatrixRealHSH
         -- * HSH Rotation Matrices
       , vecActiveRotMatrixSH
       , vecPassiveRotMatrixSH
       , vecActiveRotMatrixRealSH
       , vecPassiveRotMatrixRealSH
         -- * Test
       , testRotSymm
       , testRotSH
       --, testRotHSH
       , memoClebsh

       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U
import qualified Data.HashMap.Strict as HM
import qualified Data.Packed         as M

import           Data.Complex
import           Foreign.Storable      (Storable)
import           Data.Vector           (Vector)
import           Data.HashMap.Strict   (HashMap)
import           Data.Packed.Matrix    ((><), (@@>), cols)
import           Numeric.LinearAlgebra ((<>))
import           Numeric.Container     (Product, add, ident, ctrans, multiply, cmap)
import           System.Random

import           Texture.SH.Pyramid
import           Texture.SH.HyperSphere
import           Texture.SH.SupportFunctions
import           Texture.SphericalHarmonics
import           Texture.Orientation

import           Debug.Trace
dbg s x = trace (show s L.++ show x) x

-- =========================== Hyper Spherical Harmonics Rotation ========================

-- | Rotation matrix series for hyper spherical harmonics from 2l = 0 until 2l = 2n. More
-- info on PhD Thesis J.K. Mason eq (52).
vecRotMatrixHSH :: Int -> SO3 -> SO3 -> Vector (M.Matrix (Complex Double))
vecRotMatrixHSH n r p = let
  lmax = n `quot` 2
  us   = vecActiveRotMatrixSH  lmax r
  us'  = vecPassiveRotMatrixSH lmax p
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
    t' = M.trans $ M.mapMatrix conjugate t
    in M.mapMatrix realPart $ t' <> u <> t
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

realPartMatrixHSH :: Int -> M.Matrix (Complex Double)
realPartMatrixHSH l = let
  calcAtM (mi, mj)
    | mi == 0 && mj == 0 = (i l)
    | abs mi /= abs mj   = 0 :+ 0
    | mj < 0 && mi < 0   = (i l) * ((s / (sqrt 2)) :+ 0)
    | mj > 0 && mi < 0   = (i l) * ((1 / (sqrt 2)) :+ 0)
    | mj < 0 && mi > 0   = (i $ l-1) * ((-1 / (sqrt 2)) :+ 0)
    | otherwise          = (i $ l-1) * ((s  / (sqrt 2)) :+ 0) --  j > 0 && i > 0
    where
      i = powerComplex
      s = if even mi then 1 else (-1)
  ms = [l, l-1 .. (-l)]
  ps = [(mi, mj) | mi <- ms, mj <- ms]
  in ((2*l+1) >< (2*l+1)) $ map calcAtM ps

toRealMatrix :: M.Matrix (Complex Double) -> M.Matrix Double
toRealMatrix m = let
  l = (cols m - 1) `quot` 2
  u = realPartMatrixHSH l
  in cmap realPart $ u `multiply` m `multiply` (ctrans u)

-- ================================ Conversion to real Pyramid ===========================

-- | Transforms a complex Pyramid (N, L, MF) to its real equivalent.
fromComplexHSH :: Pyramid (N, L, MF) (Complex Double) -> Pyramid (N, L, MF) Double
fromComplexHSH pyC = let
  n       = getMaxStack pyC
  ns      = [0, 2 .. n]
  func ni = fromComplexSH (getNLSlice pyC (N ni))
  pyra    = U.concat $ map (pyramid . func) ns
  in pyC { pyramid = pyra }

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

-- ================================ Useful HSH functions =================================

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

-- | Applies a given rotation matrix to Hyper Spherical Harmonics.
applyRotMatrixHSH :: (Product a, U.Unbox a)=> Vector (M.Matrix a) -> Pyramid (N, L, MF) a -> Pyramid (N, L, MF) a
applyRotMatrixHSH rotvec = multNStack (rotvec V.!)

-- ======================================== Test ========================================

-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
-- Real SH Rotation matrix are broken. fromComplexSH doesn't work for SH but fromComplexSH_wrong works
-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

{--

rotPassiveSH + (SO3 0 0 pi) => active rotation
rotActiveSH + (SO3 0 0 (\x -> 2pi - x)) => active rotation
--}


-- | Applies both active and passive rotations on a Hyper Spherical Harmonics. The active
-- one is a rotation on the crystal reference frame and the passive is on the sample frame.
rotRHSH :: (SO3, SO3) -> Pyramid (N, L, MF) Double -> Pyramid (N, L, MF) Double
rotRHSH (ac, pa) = \hsh -> let
  n      = getMaxStack hsh
  vecmat = V.map toRealMatrix $ vecRotMatrixHSH n ac pa
  in applyRotMatrixHSH vecmat hsh

testRotSH :: IO ()
testRotSH = testSH_C (fromComplexSH_wrong . rotActiveSH (SO3 (pi/3) (pi/3) (2*pi-2*pi/5)))
--testRotSH = testSH $ rotPassiveRealSH (SO3 (-pi/2) (pi/2) (pi/2))

testRotHSH :: IO ()
testRotHSH = testHSH_C (fromComplexHSH . (rotHSH (SO3 (-pi/2) (pi/2) 0, SO3 0 0 0)))

testRotSymm2 = testHSH (rotRHSH (SO3 0 0 0, SO3 (pi/2) 0 0))

testRotSymm :: IO ()
testRotSymm = let
  symm = [ (SO3 0 0 0, SO3 (-pi/2) (pi/2) 0)
         , (SO3 (-pi/2) (pi/2) 0, SO3 0 0 0)
         , (SO3 0 0 0, SO3 0 0 0)
         ]
  rot  = symmRotMatrixHSH symm 6
  func = fromComplexHSH . applyRotMatrixHSH rot
  in testHSH_C func

testHSH_C :: (Pyramid (N,L,MF) (Complex Double) -> Pyramid (N,L,MF) Double) -> IO ()
testHSH_C func = let
  s1 = SO3 (pi/2) (pi/2) (pi/2)
  s2 = SO3 (pi/4) (pi/4) 0
  xs :: Vector (Complex Double, SO3)
  xs  = V.fromList [(10 :+ 0, s1), (10 :+ 0, s2)]
  c   = findSHCoefWeight 6 xs
  vtk = renderSO3SolidVTK (evalSH $ func c)

  rot = so3ToQuaternion $ SO3 (pi/2) (pi/2) (0)
  q1  = so3ToQuaternion s1
  q2  = so3ToQuaternion s2
  so  = U.fromList [s1, s2]
  sp  = U.fromList [ quaternionToSO3 $ q1 #<= rot
                   , quaternionToSO3 $ q2 #<= rot ]
  sa  = U.fromList [ quaternionToSO3 $ activeRot q1 rot
                   , quaternionToSO3 $ activeRot q2 rot ]
  pso = renderSO3PointsVTK so
  psp = renderSO3PointsVTK sp
  psa = renderSO3PointsVTK sa
  in do
    writeQuater "Rot-SHS-original-points" pso
    writeQuater "Rot-SHS-passive-points"  psp
    writeQuater "Rot-SHS-active-points"   psa
    writeQuater "Rot-SHS-test" vtk




testSH_C :: (Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double) -> IO ()
testSH_C func = let
  xs :: U.Vector (Complex Double, SO2)
  xs   = U.fromList [(10 :+ 0, (SO2 0 0)), (20 :+ 0, SO2 (pi/2) (pi/2))]
  pxs  = renderSO2PointsVTK $ U.map snd xs
  arxs = renderSO2PointsVTK $ U.map (arot (SO3 (pi/3) (pi/3) (2*pi/5)) . snd) xs
  prxs = renderSO2PointsVTK $ U.map (prot (SO3 (pi/3) (pi/3) (2*pi/5)) . snd) xs
  c   = findSHCoefWeight 10 (V.convert xs)
  rvtk = renderSO2VTK (evalSH  $ func c)
  ovtk = renderSO2VTK (evalSH $ fromComplexSH_wrong c)
  arot r g = cartToSO2 $ activeVecRotation  (so2ToCart g) (so3ToQuaternion r)
  prot r g = cartToSO2 $ passiveVecRotation (so2ToCart g) (so3ToQuaternion r)
  in do
    writeQuater "SH-original-points" pxs
    writeQuater "SH-active-rotated-points"  arxs
    writeQuater "SH-passive-rotated-points" prxs
    writeQuater "SH-original-ODF"    ovtk
    writeQuater "SH-rotated-ODF"     rvtk

testSH :: (Pyramid (L, MF) Double -> Pyramid (L, MF) Double) -> IO ()
testSH func = let
  xs :: U.Vector (Double, SO2)
  xs  = U.fromList [(10, SO2 0 0), (10, SO2 (pi/4) (pi/4)), (10, SO2 (pi/2) (pi/2))]
  pxs = renderSO2PointsVTK $ U.map snd xs
  rxs = renderSO2PointsVTK $ U.map (rot (SO3 (pi/2) (pi/2) 0) . snd) xs
  c   = findSHCoefWeight 10 (V.convert xs)
  rvtk = renderSO2VTK (evalSH $ func c)
  ovtk = renderSO2VTK (evalSH c)
  rot r g = cartToSO2 $ activeVecRotation (so2ToCart g) (so3ToQuaternion r)
  in do
    writeQuater "SH-original-points" pxs
    writeQuater "SH-rotated-points"  rxs
    writeQuater "SH-original-ODF"    ovtk
    writeQuater "SH-rotated-ODF"     rvtk

testSHRot :: Pyramid (L, MF) Double
testSHRot = let
  xs :: Vector (Double, SO2)
  xs = V.fromList [(10, SO2 0 0), (10, SO2 (pi/4) (pi/4)), (10, SO2 (pi/2) (pi/2))]
  in findSHCoefWeight 10 xs

testHSHRotC :: Pyramid (N, L, MF) (Complex Double)
testHSHRotC = let
  xs :: Vector (Complex Double, SO3)
  xs = V.fromList [(10 :+ 0, (SO3 0 0 0)), (10 :+ 0, SO3 (pi/2) (pi/4) 0), (10 :+ 0, SO3 (pi/2) 0 (pi/2))]
  in findSHCoefWeight 10 xs

testHSHRot :: Pyramid (N, L, MF) Double
testHSHRot = let
  xs :: Vector (Double, SO3)
  xs = V.fromList [(10, SO3 0 0 0), (10, SO3 (pi/2) (pi/4) 0), (10, SO3 (pi/2) 0 (pi/2))]
  in findSHCoefWeight 10 xs

-- ======================================== Trash ========================================

-- TODO  move to Orientation and add test case
-- | Applies an active rotation to an other rotation. It means that the second rotation is
-- applied on the reference frame instead of the rotated frame. It is used for applying
-- sample symmetry.
activeRot :: Quaternion -> Quaternion -> Quaternion
activeRot q1 q2 = let
  (q20, q2v) = splitQuaternion q2
  q2vAct     = activeVecRotation q2v (invert q1)
  in q1 #<= unsafeMergeQuaternion (q20, q2vAct)

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
