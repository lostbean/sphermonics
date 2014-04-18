{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE MultiParamTypeClasses      #-}
{-# LANGUAGE TypeFamilies               #-}

module Texture.SphericalHarmonics
       ( SH (genSHFunc)
       , evalSH
       , findSHCoef
       , findSHCoefWith
       , findSHCoefWeight
         -- * test functions (TODO remove)
       , testHSH
       , testHSH2
       , genSpread
       --, testSH
       , printSH
       , writeQuater
       , calcZ
       , calcHyperK
       , plotSHFuncFamily
       , plotHSHFuncFamily
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector         (Vector)

import           Data.Complex
import           System.Random

import           Hammer.Math.Algebra
import           Hammer.VTK

import           Texture.SH.Pyramid
import           Texture.SH.SupportFunctions
import           Texture.SH.HyperSphere

import           Debug.Trace
dbg s x = trace (s L.++ show x) x


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

{-# INLINE calcYC #-}
calcYC :: (L, MF) -> Pyramid (L, MF) Double -> Pyramid (L, MF) Double -> SO2 -> Complex Double
calcYC lmf@(_, mf) pyK pyP SO2{..}
  | mf /= 0   = kC * cis (mi * so2Phi)
  | otherwise = kC
  where
    s  = if even (unMF mf) then 1 else -1
    kC = (s * k :+ 0)
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
calcHyperK :: (N, L, MF) -> Double
calcHyperK (N n, L l, MF m) = let
  n1 = log $ fromIntegral (2 * l + 1)
  n2 = logFact (l - m)
  n3 = log $ fromIntegral (n + 1)
  n4 = logFact (n - l)
  d1 = logFact (l + m)
  d2 = logFact (n + l + 1)
  r1 = 0.5 * (n1 + n2 + n3 + n4 - d1 - d2)
  r2 = (fromIntegral l) * (log 2) + (logFact l) - (log pi)
  in exp (r1 + r2)

{-# INLINE calcZC #-}
calcZC :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, MF) Double ->
          Pyramid (N, L) Double -> SO3 -> Complex Double
calcZC nlmf@(n, l, mf) pyK pyP pyC SO3{..} = i * common * cis (mi * so3Phi)
  where
    mi = fromIntegral (unMF mf)
    li = fromIntegral (unL l) ::  Int
    s  = if even (unL l) then 1 else -1
    i  = powerComplex (unL l)
    common = let
      z1 = pyK %! nlmf
      z2 = (sin (so3Omega/2))^li
      z3 = pyC %! (n, l)
      z4 = pyP %! (l, mf)
      co = s * 0.5 * sqrt2 * z1 * z2 * z3 * z4
      in co :+ 0

{-# INLINE calcZ #-}
calcZ :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, MF) Double ->
         Pyramid (N, L) Double -> SO3 -> Double
calcZ nlmf@(n, l, mf) pyK pyP pyC SO3{..}
  | mf > 0    = s * common * cos (mi * so3Phi)
  | mf < 0    = s * common * sin (mi * so3Phi)
  | otherwise = common
  where
    mi = fromIntegral (unMF mf)
    li = fromIntegral (unL l) :: Int
    s  = if even (unMF mf) then 1 else -1
    common = let
      z1 = pyK %! nlmf
      z2 = (sin (so3Omega/2))^li
      z3 = pyC %! (n, l)
      z4 = pyP %! (l, mf)
      in z1 * z2 * z3 * z4


-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-- TESTING and IMPLEMENT MEMORIZATION

fixAngle :: Double -> Double
fixAngle x
  | x < 0     = fixAngle (x + 2*pi)
  | otherwise = x

memoLegendreFullPyramid :: Double -> Pyramid (L, MF) Double
memoLegendreFullPyramid = (mem V.!) . round . (/step) . fixAngle
  where
    mem :: Vector (Pyramid (L, MF) Double)
    mem = V.map (genAssLegenFullPyramid (L 30) . cos) ws
    ws  = V.enumFromStepN 0 step n
    n   = 360
    step = pi / (fromIntegral n)

memoGegenbauerPyramid :: Double -> Pyramid (N, L) Double
memoGegenbauerPyramid = let
  mem :: Vector (Pyramid (N, L) Double)
  mem = V.map (genGegenbauerPyramid (N 30) . cos) ws
  ws  = V.enumFromStepN 0 step n
  n   = 180
  step = pi / (fromIntegral n)
  in (mem V.!) . round . (/step) . (*0.5) . fixAngle

-- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-- =============================== Spherical Harmonic Class ==============================

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
    k  = generatePyramid calcHyperK ni
    --p  = memoLegendreFullPyramid so3Theta
    --c  = memoGegenbauerPyramid   so3Omega
    p  = genAssLegenFullPyramid (L ni) (cos so3Theta)
    c  = genGegenbauerPyramid   (N ni) (cos $ so3Omega / 2)
    in generatePyramid (\nlm -> calcZC nlm k p c r) ni

instance SH SO3 Double where
  type PyIx SO3 Double = (N, L, MF)
  genSHFunc ni = \r@(SO3{..}) -> let
    k  = generatePyramid calcHyperK ni
    p  = genAssLegenFullPyramid (L ni) (cos so3Theta)
    c  = genGegenbauerPyramid   (N ni) (cos $ so3Omega / 2)
    in generatePyramid (\nlmf -> calcZ nlmf k p c r) ni

{-# INLINE evalSH #-}
evalSH :: (SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a))=> Pyramid (PyIx pi a) a -> pi -> a
evalSH pyC = \x -> let
  ni  = getMaxStack pyC
  shf = genSHFunc ni x
  ps  = genLinSeq ni
  {-# INLINE func #-}
  func acc nlmf = acc + (pyC %! nlmf) * (shf %! nlmf)
  in L.foldl' func 0 ps

findSHCoef :: (SH pi a, U.Unbox a, Num a, PyraKey (PyIx pi a))=>
              Int -> pi -> Pyramid (PyIx pi a) a
findSHCoef = genSHFunc

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

-- ================================= Test functions ======================================

evalSingleSH :: (L, MF) -> SO2 -> Double
evalSingleSH  lmf@(l,_) = \x@(SO2{..}) -> let
  k = generatePyramid calcK (unL l)
  p = genAssLegenPyramid l (cos so2Theta)
  in calcY lmf k p x

evalSingleHSH :: (N, L, MF) -> SO3 -> Double
evalSingleHSH nlmf@(n,_,_) x@(SO3{..}) = let
  k = generatePyramid calcHyperK (unN n)
  p = genAssLegenFullPyramid (L (unN n)) (cos so3Theta)
  c = genGegenbauerPyramid n (cos $ so3Omega / 2)
  in calcZ nlmf k p c x

plotSHFuncFamily :: Int -> VTK Vec3
plotSHFuncFamily l = let
  (grid, vtk) = mkSO2 60 60
  lms         = genLinSeq l :: [(L, MF)]
  addLM :: VTK Vec3 -> (L, MF) -> VTK Vec3
  addLM acc lmf = let
    func i _ = evalSingleSH lmf (grid U.!i)
    attr = mkPointAttr (show lmf) func
    in addDataPoints acc attr
  in L.foldl addLM vtk lms

plotHSHFuncFamily :: Int -> VTK Vec3
plotHSHFuncFamily n = let
  (grid, vtk) = mkSO3 30 30 30
  lms         = genLinSeq n :: [(N, L, MF)]
  addLM :: VTK Vec3 -> (N, L, MF) -> VTK Vec3
  addLM acc nlmf = let
    func i _ = evalSingleHSH nlmf (grid U.!i)
    attr = mkPointAttr (show nlmf) func
    in addDataPoints acc attr
  in L.foldl addLM vtk lms

testZ :: SO2 -> Bool
testZ r = let
  theta = so2Theta r
  phi   = so2Phi   r
  z   = genSHFunc 7 r
  p70 = (1/32) * sqrt (15/pi)    * (429 * (cos theta)^7 - 693 * (cos theta)^5 + 315 * (cos theta)^3 - 35 * cos theta)
  p40 = (3/16) * sqrt (1/pi)     * (35 * (cos theta)^4 - 30 * (cos theta)^2 + 3)
  p30 = (1/4)  * sqrt (7/pi)     * (5 * (cos theta)^3 - 3 * cos theta)
  p32 = (1/4)  * sqrt (105/4*pi) * (cos theta * (sin theta)^2 * cos (2*phi))
  ts  = [ z %! (7, 0) - p70, z %! (4, 0) - p40, z %! (3, 0) - p30, z %! (3, 2) - p32]
  in and $ dbg "" $ map ((< 10-8) . abs) ts

testP :: Double -> Bool
testP x = let
  p = genAssLegenFullPyramid 4 x
  t1 = (p %! (3,  0)) - 0.5*(5*x^3-3*x)
  t2 = (p %! (4,  4)) - 105*(1-x^2)^2
  t3 = (p %! (4, -4)) - (105*(1-x^2)^2)/40320
  t4 = (p %! (3,  3)) - (-15*(1-x^2)**(3/2))
  t5 = (p %! (3, -3)) - (-15*(1-x^2)**(3/2))/(-720)
  in and $ dbg "" $ map ((< 10-8) . abs) [t1, t2, t3, t4, t5]

testSH :: (Pyramid (L, MF) Double -> Pyramid (L, MF) Double) -> IO ()
testSH func = let
  xs :: Vector (Double, SO2)
  xs  = V.fromList [(10, SO2 0 (pi/4)), (10, SO2 (pi/4) (pi/2)), (10, SO2 (pi/2) (3*pi/4))]
  --xs  = V.fromList [(10, (0, pi/4)), (10, (pi/4, pi/2)), (10, (pi/2, 3*pi/4))]
  c   = findSHCoefWeight 10 xs
  vtk = renderSO2VTK (evalSH $ func c)
  in writeQuater "ODF-SH-test-" vtk

genSpread :: SO3 -> IO (Vector (Double, SO3))
genSpread SO3{..} = let
  s = 10*pi/180
  func = do
    o <- randomRIO (so3Omega - s, so3Omega + s)
    t <- randomRIO (so3Theta - s, so3Theta + s)
    p <- randomRIO (so3Phi   - s, so3Phi   + s)
    return (1, SO3 o t p)
  in V.replicateM 50 func

testHSH :: (Pyramid (N,L,MF) Double -> Pyramid (N,L,MF) Double) -> IO ()
testHSH func = let
  xs :: Vector (Double, SO3)
  xs  = V.fromList [(10, SO3 (pi/2) (pi/2) (pi/2)), (10, SO3 pi (pi/4) 0)]
  c   = findSHCoefWeight 7 xs
  vtk = renderSO3SolidVTK (evalSH $ func c)
  in writeUniVTKfile ("/home/edgar/Desktop/SHS-test.vtu") False vtk

testHSH2 :: Int -> Vector (Double, SO3) -> IO ()
testHSH2 n xs = let
  c   = findSHCoefWeight n xs
  ps  = renderSO3PointsVTK $ V.convert $ V.map snd xs
  vtk = renderSO3SolidVTK (evalSH c)
  in do
    writeUniVTKfile ("/home/edgar/Desktop/SHS-test-points.vtu") False ps
    writeUniVTKfile ("/home/edgar/Desktop/SHS-test.vtu") False vtk

printSH :: Pyramid (L, MF) (Double) -> IO ()
printSH c = let
  vtk = renderSO2VTK (evalSH c)
  in writeUniVTKfile ("/home/edgar/Desktop/SH-test.vtu") False vtk

printHSH :: Pyramid (N, L, MF) (Double) -> IO ()
printHSH c = let
  vtk = renderSO3SolidVTK (evalSH c)
  in writeUniVTKfile ("/home/edgar/Desktop/HSH-test.vtu") False vtk

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") False
