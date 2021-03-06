{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE MultiParamTypeClasses      #-}
{-# LANGUAGE TypeFamilies               #-}

module Texture.SH.Harmonics
       ( SH (genSHFunc)
       , evalSH
       , findSHCoef
       , findSHCoefWith
       , findSHCoefWeight
         -- * Plotting
       , plotSH
       , plotHSH
       , plotSH_C
       , plotHSH_C
       , plotSHPoints
       , plotHSHPoints
       , plotSHFuncFamily
       , plotHSHFuncFamily
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector         (Vector)
import           System.IO.Unsafe    (unsafePerformIO)
import           Control.Exception   (try, SomeException)
import           Data.Binary

import           Data.Complex

import           Hammer.Math.Algebra
import           Hammer.VTK

import           Texture.HyperSphere
import           Texture.Orientation

import           Texture.SH.Pyramid
import           Texture.SH.SupportFunctions

--import           Debug.Trace
--dbg s x = trace (s L.++ show x) x

sqrt2 :: Double
sqrt2 = sqrt 2

-- ================================= Spherical Harmonics =================================

{-# INLINE calcK #-}
calcK :: (L, M) -> Double
calcK (l, m) = calcKFull (l, m2mf m)

calcKFull :: (L, MF) -> Double
calcKFull (L l, MF mf) = let
  a = log $ fromIntegral $ 2 * l + 1
  b = logFact (l - mf)
  c = logFact (l + mf)
  in sqrt $ exp (a + b - log (4*pi) - c)

-- TODO use conjugate and reduce Pyramid size (L, MF) -> (L, M) for pyK and pyP
{-# INLINE calcYC #-}
calcYC :: (L, MF) -> Pyramid (L, MF) Double -> Pyramid (L, MF) Double -> SO2 -> Complex Double
calcYC lmf@(_, mf) pyK pyP SO2{..}
  | mf /= 0   = kC * cis (mi * so2Phi)
  | otherwise = kC
  where
    -- Condon-Shortley phase already included in P (Legendre)
    kC = (k :+ 0)
    k  = pyK %! lmf * pyP %! lmf
    mi = fromIntegral (unMF mf)

{-# INLINE calcY #-}
calcY :: (L, MF) -> Pyramid (L, M) Double -> Pyramid (L, M) Double -> SO2 -> Double
calcY (l, mf) pyK pyP SO2{..}
  | mf > 0    = s * sqrt2 * k * cos ( mi * so2Phi)
  | mf < 0    = s * sqrt2 * k * sin (-mi * so2Phi)
  | otherwise = k
  where
    s  = if even (unMF mf) then 1 else -1
    k  = (pyK %! lm) * pyP %! lm
    lm = (l, mf2m mf)
    mi = fromIntegral (unMF mf)

-- =============================== HyperSpherical Harmonics ==============================

{-# INLINE calcHyperK #-}
calcHyperK :: (N, L, M) -> Double
calcHyperK (n, l, m) = calcHyperKFull (n, l, m2mf m)

-- | Constant used in the hyperspherical function
calcHyperKFull :: (N, L, MF) -> Double
calcHyperKFull (N n, L l, MF m) = let
  n1 = log $ fromIntegral (2 * l + 1)
  n2 = logFact (l - m)
  n3 = log $ fromIntegral (n + 1)
  n4 = logFact (n - l)
  d1 = logFact (l + m)
  d2 = logFact (n + l + 1)
  r1 = 0.5 * (n1 + n2 + n3 + n4 - d1 - d2)
  r2 = (fromIntegral l) * (log 2) + (logFact l) - (log pi)
  in exp (r1 + r2)

-- | Complex hyperspherical harmonic function for pre-calculated constant, Legendre and
-- Gegenbauer.
calcZC :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, MF) Double ->
          Pyramid (N, L) Double -> SO3 -> Complex Double
calcZC nlmf@(n, l, mf) pyK pyP pyC SO3{..} = i * common * cis (mi * so3Phi)
  where
    mi = fromIntegral (unMF mf)
    li = fromIntegral (unL l) ::  Int
    i  = negativeImaginaryPower (unL l)
    common = let
      z1 = pyK %! nlmf
      z2 = (sin (so3Omega/2))^li
      z3 = pyC %! (n, l)
      z4 = pyP %! (l, mf)
      co = 0.5 * sqrt2 * z1 * z2 * z3 * z4
      in co :+ 0

-- | Real hyperspherical harmonic function for pre-calculated constant, Legendre and
-- Gegenbauer.
calcZ :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, MF) Double ->
         Pyramid (N, L) Double -> SO3 -> Double
calcZ (n, l, mf) pyK pyP pyC x@SO3{..}
  | mf > 0    = realPart $ (ip l)     * kr * (om * zP + zN)
  | mf < 0    = realPart $ (ip $ l-1) * kr * (om * zP - zN)
  | otherwise = realPart $ (ip l)     * zP
  where
    zP  = calcZC (n, l, abs mf) pyK pyP pyC x
    zN  = olm * conjugate zP
    ip  = imaginaryPower . unL
    kr  = (1 / sqrt2) :+ 0
    om  = if even (unMF mf) then 1 else -1
    olm = if even (abs (unL l) + abs (unMF mf)) then 1 else -1

-- | This implementation has a shifted phase. See the PhD Thesis "Analysis of
-- crystallographic texture information by the hyperspherical harmonic expansion", pages
-- 137.
calcZ_bad :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, MF) Double ->
             Pyramid (N, L) Double -> SO3 -> Double
calcZ_bad (n, l, mf) pyK pyP pyC SO3{..}
  | mf > 0    = s * common * cos (mi * so3Phi)
  | mf < 0    = s * common * sin (mi * so3Phi)
  | otherwise = common
  where
    mi = fromIntegral (unMF mf)
    li = fromIntegral (unL l) :: Int
    s  = if even (unMF mf) then 1 else -1
    common = let
      z1 = pyK %! (n, l, mf)
      z2 = (sin (so3Omega/2))^li
      z3 = pyC %! (n, l)
      z4 = pyP %! (l, mf)
      in z1 * z2 * z3 * z4

-- ====================================== Memorization ===================================

fixAngle :: Double -> Double
fixAngle x
  | x < 0     = fixAngle (x + 2*pi)
  | x > 2*pi  = fixAngle (x - 2*pi)
  | otherwise = x


filetable :: (Binary a)=> FilePath -> Maybe (Vector a)
filetable file = unsafePerformIO $ do
  t <- try (decodeFile file)
  case t of
   Right d -> return $ Just d
   Left  e -> do
     putStr $ "[Harmonics] Couldn't load the data table on " ++ file
     putStrLn $ ". Reason: " ++ show (e :: SomeException)
     return Nothing

memoLegendreFullPyramid :: Double -> Pyramid (L, MF) Double
memoLegendreFullPyramid = (mem V.!) . round . (/step) . fixAngle
  where
    mem :: Vector (Pyramid (L, MF) Double)
    mem = maybe
          (V.map (genAssLegenFullPyramid (L 30) . cos) ws)
          id
          (filetable "legendre_coef.data")
    ws  = V.enumFromStepN 0 step n
    n   = 360
    step = (2 * pi) / (fromIntegral n)

memoGegenbauerPyramid :: Double -> Pyramid (N, L) Double
memoGegenbauerPyramid = (mem V.!) . round . (/step) . (* 0.5) . fixAngle
  where
    mem :: Vector (Pyramid (N, L) Double)
    mem = maybe
          (V.map (genGegenbauerPyramid (N 30) . cos) ws)
          id
          (filetable "gegenbauer_coef.data")
    ws  = V.enumFromStepN 0 step n
    n   = 180
    step = pi / (fromIntegral n)

-- =============================== Spherical Harmonic Class ==============================

-- | Class of spherical harmonics coefficient generator
class SH a b where
  type PyIx a b :: *
  genSHFunc :: Int -> a -> Pyramid (PyIx a b) b

instance SH SO2 (Complex Double) where
  type PyIx SO2 (Complex Double) = (L, MF)
  genSHFunc li = \r@(SO2{..}) -> let
    k = generatePyramid calcKFull li
    p = genAssLegenFullPyramid (L li) (cos so2Theta)
    in generatePyramid (\lmf -> calcYC lmf k p r) li

instance SH SO2 Double where
  type PyIx SO2 Double = (L, MF)
  genSHFunc li = \r@(SO2{..}) -> let
    k = generatePyramid calcK li
    p = genAssLegenPyramid (L li) (cos so2Theta)
    in generatePyramid (\lmf -> calcY lmf k p r) li

instance SH SO3 (Complex Double) where
  type PyIx SO3 (Complex Double) = (N, L, MF)
  genSHFunc ni = \r@(SO3{..}) -> let
    k  = generatePyramid calcHyperKFull ni
    p  = memoLegendreFullPyramid so3Theta
    c  = memoGegenbauerPyramid   so3Omega
    --p  = genAssLegenFullPyramid (L ni) (cos so3Theta)
    --c  = genGegenbauerPyramid   (N ni) (cos $ so3Omega / 2)
    in generatePyramid (\nlm -> calcZC nlm k p c r) ni

instance SH SO3 Double where
  type PyIx SO3 Double = (N, L, MF)
  genSHFunc ni = \r@(SO3{..}) -> let
    k  = generatePyramid calcHyperKFull ni
    p  = memoLegendreFullPyramid so3Theta
    c  = memoGegenbauerPyramid   so3Omega
    --p  = genAssLegenFullPyramid (L ni) (cos so3Theta)
    --c  = genGegenbauerPyramid   (N ni) (cos $ so3Omega / 2)
    in generatePyramid (\nlmf -> calcZ nlmf k p c r) ni

-- | Evaluate the spherical harmonic expansion represented by a given set of coefficients
-- at a given position.
{-# INLINE evalSH #-}
evalSH :: (SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a))=> Pyramid (PyIx pi a) a -> pi -> a
evalSH pyC = \x -> let
  ni  = getMaxStack pyC
  shf = genSHFunc ni x
  ps  = genLinSeq ni
  {-# INLINE func #-}
  func acc nlmf = acc + (pyC %! nlmf) * (shf %! nlmf)
  in L.foldl' func 0 ps

-- | Generate the spherical harmonics coefficient for a single position.
findSHCoef :: (SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a))=>
              Int -> pi -> Pyramid (PyIx pi a) a
findSHCoef = genSHFunc

-- | Generate the spherical harmonics coefficient for multiple positions with weights.
findSHCoefWeight :: ( SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a)
                    , Eq (PyIx pi a), Fractional a)=>
                    Int -> Vector (a, pi) -> Pyramid (PyIx pi a) a
findSHCoefWeight n xs = let
  shf = genSHFunc n
  np  = fromIntegral $ V.length xs
  x0  = V.head xs
  xt  = V.tail xs
  foo (w, x) = mapPyramid (w *) $ shf x
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/np) total

-- | Generate the spherical harmonics coefficient for multiple positions with custom
-- weight function.
findSHCoefWith :: ( SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a)
                  , Eq (PyIx pi a), Fractional a)=>
                  Int -> (pi -> a) -> Vector pi -> Pyramid (PyIx pi a) a
findSHCoefWith n func xs = let
  shf   = genSHFunc n
  np    = fromIntegral $ V.length xs
  x0    = V.head xs
  xt    = V.tail xs
  foo x = mapPyramid ((func x) *) $ shf x
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/np) total

-- =============================== Evaluating Functions ==================================

-- | Evaluate a given SH base function at a given position.
evalSingleSH_C :: (L, MF) -> SO2 -> Complex Double
evalSingleSH_C lmf@(l,_) = \x@(SO2{..}) -> let
  k = generatePyramid calcKFull (unL l)
  p = genAssLegenFullPyramid l (cos so2Theta)
  in calcYC lmf k p x

-- | Evaluate a given SH base function at a given position. Uses 'calcZ'
evalSingleSH :: (L, MF) -> SO2 -> Double
evalSingleSH  lmf@(l,_) = \x@(SO2{..}) -> let
  k = generatePyramid calcK (unL l)
  p = genAssLegenPyramid l (cos so2Theta)
  in calcY lmf k p x

-- | Evaluate a given SH base function at a given position. Uses 'calcZC'
evalSingleSH' :: (L, MF) -> SO2 -> Double
evalSingleSH' lmf@(l, mf) x@(SO2{..})
  | mf > 0    = s * realPart (func (l, mf) + conjugate (func (l, mf)))
  | mf < 0    = s * realPart (ni * (func (l, -mf) - (conjugate $ func (l, -mf))))
  | otherwise = realPart $ func (l, mf)
  where
    s  = if even (unMF mf) then sr else -sr
    sr = (1 / (sqrt 2))
    ni = (0 :+ (-1))
    k  = generatePyramid calcKFull (unL l)
    p  = genAssLegenFullPyramid l (cos so2Theta)
    func pos = calcYC pos k p x

-- | Evaluate a given SH base function at a given position. Uses 'calcZC'
evalSingleSH'' :: (L, MF) -> SO2 -> Double
evalSingleSH'' lmf@(l, mf) x@(SO2{..})
  | mf > 0    = realPart $ kr * (s * func (l, mf)  + (func (l, -mf)))
  | mf < 0    = realPart $ ki * (s * func (l, -mf) - (func (l,  mf)))
  | otherwise = realPart $ func (l, mf)
  where
    kr = (1 / sqrt2) :+ 0
    ki = 0 :+ (-1 / sqrt2)
    s  = if even (unMF mf) then 1 else -1
    k  = generatePyramid calcKFull (unL l)
    p  = genAssLegenFullPyramid l (cos so2Theta)
    func pos = calcYC pos k p x

-- | Evaluate a given HSH base function at a given position. Uses 'calcZ'
evalSingleHSH :: (N, L, MF) -> SO3 -> Double
evalSingleHSH nlmf@(n,_,_) x@(SO3{..}) = let
  k = generatePyramid calcHyperKFull (unN n)
  p = genAssLegenFullPyramid (L (unN n)) (cos so3Theta)
  c = genGegenbauerPyramid n (cos $ so3Omega / 2)
  in calcZ nlmf k p c x

-- | Evaluate a given HSH base function at a given position. Uses 'calcZC'
evalSingleHSH' :: (N, L, MF) -> SO3 -> Double
evalSingleHSH' (n, l, mf) x@(SO3{..})
  | mf > 0    = realPart $ (ip l)     * kr * (s * func (n, l, mf)  + (func (n, l, -mf)))
  | mf < 0    = realPart $ (ip $ l-1) * kr * (s * func (n, l, -mf) - (func (n, l,  mf)))
  | otherwise = realPart $ (ip l)     * func (n, l, mf)
  where
    ip = imaginaryPower . unL
    kr = (1 / sqrt2) :+ 0
    s  = if even (unMF mf) then 1 else -1
    k = generatePyramid calcHyperKFull (unN n)
    p = genAssLegenFullPyramid (L (unN n)) (cos so3Theta)
    c = genGegenbauerPyramid n (cos $ so3Omega / 2)
    func pos = calcZC pos k p c x

-- =============================== Plotting Functions ====================================

-- | Calculates and plots SH given a list of positions and a function to convert from
-- Complex SH to Real SH.
plotSH_C :: Int
         -> [SO2]
         -> (Pyramid (L, MF) (Complex Double) -> Pyramid (L, MF) Double)
         -> VTK Vec3
plotSH_C l ss func = let
  xs :: Vector (Complex Double, SO2)
  xs  = V.fromList $ map (\s -> (1 :+ 0, s)) ss
  c   = findSHCoefWeight l xs
  in renderSO2VTK (evalSH (func c))

-- | Calculates and plots SH given a list of positions and a transformation function (e.g.
-- rotation) to convert from Complex SH to Real SH.
plotSH :: Int
       -> [SO2]
       -> (Pyramid (L, MF) Double -> Pyramid (L, MF) Double)
       -> VTK Vec3
plotSH l ss func = let
  xs :: Vector (Double, SO2)
  xs  = V.fromList $ map (\s -> (1, s)) ss
  c   = findSHCoefWeight l xs
  in renderSO2VTK (evalSH $ func c)

-- | Calculates and plots HSH given a list of positions and a function to convert from
-- Complex SH to Real SH.
plotHSH_C :: Int
          -> [SO3]
          -> (Pyramid (N,L,MF) (Complex Double) -> Pyramid (N,L,MF) Double)
          -> VTK Vec3
plotHSH_C n ss func = let
  xs :: Vector (Complex Double, SO3)
  xs  = V.fromList $ map (\s -> (1 :+ 0, s)) ss
  c   = findSHCoefWeight n xs
  in renderSO3SolidVTK (evalSH (func c))

-- | Calculates and plots HSH given a list of positions and a transformation function (e.g.
-- rotation) to convert from Complex SH to Real SH.
plotHSH :: Int
        -> [SO3]
        -> (Pyramid (N,L,MF) Double -> Pyramid (N,L,MF) Double)
        -> VTK Vec3
plotHSH n ss func = let
  xs :: Vector (Double, SO3)
  xs  = V.fromList $ map (\s -> (1, s)) ss
  c   = findSHCoefWeight n xs
  in renderSO3SolidVTK (evalSH $ func c)

-- ========================= Plotting Functions for testing ==============================

-- | Plots all base functions up two a give L.
plotSHFuncFamily :: Int -> IO ()
plotSHFuncFamily l = let
  (grid, vtk) = mkSO2 60 60
  lms         = genLinSeq l :: [(L, MF)]
  addLM :: VTK Vec3 -> (L, MF) -> VTK Vec3
  addLM acc lmf = let
    func = evalSingleSH lmf . (grid U.!)
    attr = mkPointAttr (show lmf) func
    in addPointAttr acc attr
  out = L.foldl addLM vtk lms
  in writeQuater ("SH_Family-L=" ++ show l) out

-- | Plots all base functions up two a give N.
plotHSHFuncFamily :: Int -> IO ()
plotHSHFuncFamily n = let
  (grid, vtk) = mkSO3 30 30 30
  lms         = genLinSeq n :: [(N, L, MF)]
  addLM :: VTK Vec3 -> (N, L, MF) -> VTK Vec3
  addLM acc nlmf = let
    func = evalSingleHSH nlmf . (grid U.!)
    attr = mkPointAttr (show nlmf) func
    in addPointAttr acc attr
  out = L.foldl addLM vtk lms
  in writeQuater ("HSH_Family-N=" ++ show n) out

-- | Plots lists of positions and their active and passive rotations in the SO2 domain.
plotSHPoints :: [SO2] -> [SO3] -> [SO3] -> IO ()
plotSHPoints ss as ps = let
  ss'   = U.fromList ss
  as'   = U.map so3ToQuaternion $ U.fromList as
  ps'   = U.map so3ToQuaternion $ U.fromList ps
  qplot = renderSO2PointsVTK ss'
  aplot = renderSO2PointsVTK $ U.concatMap (\r -> U.map (arot r) ss') as'
  pplot = renderSO2PointsVTK $ U.concatMap (\r -> U.map (prot r) ss') ps'
  arot r g = cartToSO2 $ activeVecRotation  (so2ToCart g) r
  prot r g = cartToSO2 $ passiveVecRotation (so2ToCart g) r
  in do
    writeQuater "SH-original-points" qplot
    writeQuater "SH-active-points"   aplot
    writeQuater "SH-passive-points"  pplot

-- | Plots lists of positions and their active and passive rotations in the SO3 domain.
plotHSHPoints :: [SO3] -> [SO3] -> [SO3] -> IO ()
plotHSHPoints ss as ps = let
  ss'   = U.map so3ToQuaternion $ U.fromList ss
  as'   = U.map so3ToQuaternion $ U.fromList as
  ps'   = U.map so3ToQuaternion $ U.fromList ps
  qplot = renderSO3PointsVTK $ U.fromList ss
  aplot = renderSO3PointsVTK $ U.concatMap (\r -> U.map (arot r) ss') as'
  pplot = renderSO3PointsVTK $ U.concatMap (\r -> U.map (prot r) ss') ps'
  arot r g = quaternionToSO3 $ g #<= r
  prot r g = quaternionToSO3 $ activeRot g r
  in do
    writeQuater "HSH-original-points" qplot
    writeQuater "HSH-active-points"   aplot
    writeQuater "HSH-passive-points"  pplot

-- TODO  move to Orientation and add test case
-- | Applies an active rotation to an other rotation. It means that the second rotation is
-- applied on the reference frame instead of the rotated frame. It is used for applying
-- sample symmetry.
activeRot :: Quaternion -> Quaternion -> Quaternion
activeRot q1 q2 = let
  (q20, q2v) = splitQuaternion q2
  q2vAct     = activeVecRotation q2v (invert q1)
  in q1 #<= unsafeMergeQuaternion (q20, q2vAct)

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") False

-- ================================= Test functions ======================================

testZ :: SO2 -> Bool
testZ r = let
  theta = so2Theta r
  phi   = so2Phi   r
  z   = genSHFunc 7 r
  p70 = (1/32) * sqrt (15/pi)    * (429 * (cos theta)^(7::Int) -
                                    693 * (cos theta)^(5::Int) +
                                    315 * (cos theta)^(3::Int) -
                                    35  *  cos theta)
  p40 = (3/16) * sqrt (1/pi)     * (35 * (cos theta)^(4::Int)-
                                    30 * (cos theta)^(2::Int)+ 3)
  p30 = (1/4)  * sqrt (7/pi)     * (5 *  (cos theta)^(3::Int)-
                                    3 *   cos theta)
  p32 = (1/4)  * sqrt (105/4*pi) * (cos theta * (sin theta)^(2::Int) * cos (2*phi))
  ts  = [ z %! (7, 0) - p70, z %! (4, 0) - p40, z %! (3, 0) - p30, z %! (3, 2) - p32]
  in and $ map ((< 10-8) . abs) ts

testP :: Double -> Bool
testP x = let
  p = genAssLegenFullPyramid 4 x
  t1 = (p %! (3,  0)) - 0.5*(5*x^(3::Int)-3*x)
  t2 = (p %! (4,  4)) - 105*(1-x^(2::Int))^(2::Int)
  t3 = (p %! (4, -4)) - (105*(1-x^(2::Int))^(2::Int))/40320
  t4 = (p %! (3,  3)) - (-15*(1-x^(2::Int))**(3/2))
  t5 = (p %! (3, -3)) - (-15*(1-x^(2::Int))**(3/2))/(-720)
  in and $ map ((< 10-8) . abs) [t1, t2, t3, t4, t5]

testSHGrid ::  U.Vector (Bool, Bool)
testSHGrid = let
  (grid, _) = mkSO2 30 30
  in U.map testSH grid

testSH :: SO2 -> (Bool, Bool)
testSH x@SO2{..} = let
  y00  = (0.5 / sqrt pi) :+ 0
  y1n1 = ((0.5  * sqrt (3/(2*pi))   * sin so2Theta) :+ 0) * cis (-so2Phi)
  y10  = ((0.5  * sqrt (3/pi)       * cos so2Theta) :+ 0)
  y11  = ((-0.5 * sqrt (3/(2*pi))   * sin so2Theta) :+ 0) * cis so2Phi
  y2n2 = ((0.25 * sqrt (15/(2*pi))  * (sin so2Theta)^(2::Int)) :+ 0) * cis (-2 * so2Phi)
  y2n1 = ((0.5  * sqrt (15/(2*pi))  * sin so2Theta* cos so2Theta) :+ 0) * cis (-so2Phi)
  y20  = ((0.25 * sqrt (5/pi)       * (3 * (cos so2Theta)^(2::Int) - 1)) :+ 0)
  y21  = ((-0.5 * sqrt (15/(2*pi))  * sin so2Theta* cos so2Theta) :+ 0) * cis (so2Phi)
  y22  = ((0.25 * sqrt (15/(2*pi))  * (sin so2Theta)^(2::Int)) :+ 0) * cis (2 * so2Phi)
  t1 = evalSingleSH_C (0, 0) x - y00
  t2 = evalSingleSH_C (1,-1) x - y1n1
  t3 = evalSingleSH_C (1, 0) x - y10
  t4 = evalSingleSH_C (1, 1) x - y11
  t5 = evalSingleSH_C (2,-2) x - y2n2
  t6 = evalSingleSH_C (2,-1) x - y2n1
  t7 = evalSingleSH_C (2, 0) x - y20
  t8 = evalSingleSH_C (2, 1) x - y21
  t9 = evalSingleSH_C (2, 2) x - y22

  s = sqrt 0.5
  r00  = realPart y00
  r1n1 = realPart $ (0 :+ s) * (y1n1 + y11)
  r10  = realPart $ y10
  r11  = realPart $ (s :+ 0) * (y1n1 - y11)
  r2n2 = realPart $ (0 :+ s) * (y2n2 - y22)
  r2n1 = realPart $ (0 :+ s) * (y2n1 + y21)
  r20  = realPart $ y20
  r21  = realPart $ (s :+ 0) * (y2n1 - y21)
  r22  = realPart $ (s :+ 0) * (y2n2 + y22)
  tr1 = evalSingleSH (0, 0) x - r00
  tr2 = evalSingleSH (1,-1) x - r1n1
  tr3 = evalSingleSH (1, 0) x - r10
  tr4 = evalSingleSH (1, 1) x - r11
  tr5 = evalSingleSH (2,-2) x - r2n2
  tr6 = evalSingleSH (2,-1) x - r2n1
  tr7 = evalSingleSH (2, 0) x - r20
  tr8 = evalSingleSH (2, 1) x - r21
  tr9 = evalSingleSH (2, 2) x - r22

  testC = and $ map ((< 10e-8) . magnitude) [t1, t2, t3, t4, t5, t6, t7, t8, t9]
  testR = and $ map ((< 10e-8) . abs) [tr1, tr2, tr3, tr4, tr5, tr6, tr7, tr8, tr9]
  in (testC, testR)
