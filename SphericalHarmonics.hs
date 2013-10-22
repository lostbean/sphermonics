{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Hammer.Texture.SphericalHarmonics
       ( evalRealSH
       , findRealSHCoef
       , findRealSHCoefWith
       , findRealSHCoefWeight
         -- * test functions (TODO remove)
       , testHSHSingle
       , testHSH
       , writeQuater
       , calcZ
       , calcHyperK
       , genYPyramid
       , genZPyramid
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector         (Vector)

import           Hammer.Math.Algebra
import           Hammer.Render.VTK.VTKRender
import           Hammer.Texture.Orientation
import           Hammer.Texture.SH.Pyramid
import           Hammer.Texture.SH.SupportFunctions

import           Debug.Trace
dbg s x = trace (s L.++ show x) x

-- ================================= Spherical Harmonics =================================

sqrt2 :: Double
sqrt2 = sqrt 2

{-# INLINE calcK #-}
calcK :: (L, MF) -> Double
calcK (L l, MF m) = let
  a = log $ fromIntegral $ 2 * l + 1
  b = logFact (l - abs m)
  c = logFact (l + abs m)
  in sqrt $ exp (a + b - log (4*pi) - c)

{-# INLINE calcY #-}
calcY :: (L, MF) -> Pyramid (L, MF) Double -> Pyramid (L, M) Double -> (Double, Double) -> Double
calcY lmf@(l, mf) pyK pyP (_, phi)
  | mf > 0    = sqrt2 * k * cos ( mi * phi)
  | mf < 0    = sqrt2 * k * sin (-mi * phi)
  | otherwise = k
  where
    k  = (pyK %! lmf) * pyP %! lm
    lm = (l, mf2m mf)
    mi = fromIntegral (unMF mf)

genYPyramid :: Int -> (Double, Double) -> Pyramid (L, MF) Double
genYPyramid li = \r@(theta, _) -> let
  k = generatePyramid calcK li
  p = genAssLegenPyramid (L li) (cos theta)
  in generatePyramid (\lmf -> calcY lmf k p r) li

{-# INLINE evalRealSH #-}
evalRealSH :: Pyramid (L, MF) Double -> (Double, Double) -> Double
evalRealSH pyC = let
  li = maxStack pyC
  k  = generatePyramid calcK li
  ps = genLinSeq li
  in \x@(theta, _) -> let
    p = genAssLegenPyramid (L li) (cos theta)
    func acc lmf = acc + (pyC %! lmf) * (calcY lmf k p x)
    in L.foldl' func 0 ps

findRealSHCoef :: Int -> (Double, Double) -> Pyramid (L, MF) Double
findRealSHCoef l x@(theta, _) = let
  p = genAssLegenPyramid (L l) (cos theta)
  k = generatePyramid calcK l
  func lmf = calcY lmf k p x
  in generatePyramid func l

findRealSHCoefWeight :: Int -> Vector (Double, (Double, Double)) -> Pyramid (L, MF) Double
findRealSHCoefWeight l xs = let
  k = generatePyramid calcK l
  n = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo (w, x@(theta, _)) = let
    p = genAssLegenPyramid (L l) (cos theta)
    in generatePyramid (\lmn -> w * calcY lmn k p x) l
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/n) total

findRealSHCoefWith :: Int -> ((Double, Double) -> Double) ->
                      Vector (Double, Double) -> Pyramid (L, MF) Double
findRealSHCoefWith l func xs = let
  k = generatePyramid calcK l
  n = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo x@(theta, _) = let
    p = genAssLegenPyramid (L l) (cos theta)
    f = func x
    in generatePyramid (\lmf -> f * calcY lmf k p x) l
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/n) total

-- =============================== HyperSpherical Harmonics ==============================

{-# INLINE calcHyperK #-}
calcHyperK :: (N, L, MF) -> Double
calcHyperK (N n, L l, MF m) = let
  a = log $ fromIntegral $ 2 * l + 1
  b = logFact (l - abs m)
  c = logFact (l + abs m)
  d = logFact (n - l)
  e = logFact (n + l + 1)
  r = 0.5 * (a + b + log (fromIntegral $ n+1) + d - c - e)
  in exp $ (fromIntegral l) * (log 2) + (logFact l) + r - (log pi)

{-# INLINE calcZ #-}
calcZ :: (N, L, MF) -> Pyramid (N, L, MF) Double -> Pyramid (L, M) Double ->
         Pyramid (N, L) Double -> (Double, Double, Double) -> Double
calcZ nlmf@(n, l, mf) pyK pyP pyC (omega, _ , phi)
  | mf > 0    = sign * common * cos (mi * phi)
  | mf < 0    = sign * common * sin (mi * phi)
  | otherwise = common
  where
    sign   = if even (unMF mf) then 1 else -1
    common = let
      z1 = testNan "z1" $ pyK %! nlmf
      z2 = testNan "z2" $ (sin (omega/2))^li
      z3 = testNan "z3" $ pyC %! (n, l)
      z4 = testNan (show nlmf) $ pyP %! (l, mf2m mf)
      in z1 * z2 * z3 * z4
    mi     = fromIntegral (unMF mf)
    li     = fromIntegral (unL l) ::  Int

genZPyramid :: Int -> (Double, Double, Double) -> Pyramid (N, L, MF) Double
genZPyramid ni = \r@(omega, theta, _) -> let
  k  = generatePyramid calcHyperK ni
  p  = genAssLegenPyramid   (L ni) (cos theta)
  c  = genGegenbauerPyramid (N ni) (cos $ omega / 2)
  in generatePyramid (\nlmf -> calcZ nlmf k p c r) ni

testNan s x
  | isNaN x   = error $ s ++ " - " ++ show x
  | otherwise = x

{-# INLINE evalRealHSH #-}
evalRealHSH :: Pyramid (N, L, MF) Double -> (Double, Double, Double) -> Double
evalRealHSH pyC  = let
  ni = maxStack pyC
  k  = generatePyramid calcHyperK ni
  ps = genLinSeq ni
  in \x@(omega, theta, _) -> let
    {-# INLINE func #-}
    func acc nlmf = acc + (pyC %! nlmf) * (calcZ nlmf k p c x)
    p  = genAssLegenPyramid   (L ni) (cos theta)
    c  = genGegenbauerPyramid (N ni) (cos $ omega / 2)
    in L.foldl' func 0 ps

findRealHSHCoef :: Int -> (Double, Double, Double) -> Pyramid (N, L, MF) Double
findRealHSHCoef n x@(omega, theta, _) = let
  p = genAssLegenPyramid (L n) (cos theta)
  k = generatePyramid calcHyperK n
  c = genGegenbauerPyramid (N n) (cos $ omega / 2)
  func nlmf = calcZ nlmf k p c x
  in generatePyramid func n

findRealHSHCoefWeight :: Int -> Vector (Double, (Double, Double, Double)) ->
                         Pyramid (N, L, MF) Double
findRealHSHCoefWeight n xs = let
  k  = generatePyramid calcHyperK n
  np = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo (w, x@(omega, theta, _)) = let
    p = genAssLegenPyramid   (L n) (cos theta)
    c = genGegenbauerPyramid (N n) (cos $ omega / 2)
    in generatePyramid (\nlmf -> w * calcZ nlmf k p c x) n
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/np) total

findRealHSHCoefWith :: Int -> ((Double, Double, Double) -> Double) ->
                       Vector (Double, Double, Double) -> Pyramid (N, L, MF) Double
findRealHSHCoefWith n func xs = let
  k  = generatePyramid calcHyperK n
  np = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo x@(omega, theta, _) = let
    p = genAssLegenPyramid   (L n) (cos theta)
    c = genGegenbauerPyramid (N n) (cos $ omega / 2)
    f = func x
    in generatePyramid (\nlmf -> f * calcZ nlmf k p c x) n
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/np) total

-- ================================= Test functions ======================================

testSHSingle :: (L, MF) -> VTK Vec3
testSHSingle lmf = renderSphereVTK (evalSingleSH lmf)

evalSingleSH :: (L, MF) -> (Double, Double) -> Double
evalSingleSH  lmf@(l,_) x@(theta, _) = let
  k = generatePyramid calcK (unL l)
  p = genAssLegenPyramid l (cos theta)
  in calcY lmf k p x

testHSHSingle :: (N, L, MF) -> VTK Vec3
testHSHSingle nlmf = renderHyperSphereVTK (evalSingleHSH nlmf)

evalSingleHSH :: (N, L, MF) -> (Double, Double, Double) -> Double
evalSingleHSH nlmf@(n, l,_) x@(omega, theta, _) = let
  k = generatePyramid calcHyperK (unN n)
  p = genAssLegenPyramid (L (unN n)) (cos theta)
  c = genGegenbauerPyramid n (cos $ omega / 2)
  in calcZ nlmf k p c x

p70 (theta, phi) = (1/32) * sqrt (15/pi)    * (429 * (cos theta)^7 - 693 * (cos theta)^5 + 315 * (cos theta)^3 - 35 * cos theta)
p40 (theta, phi) = (3/16) * sqrt (1/pi)     * (35 * (cos theta)^4 - 30 * (cos theta)^2 + 3)
p30 (theta, phi) = (1/4)  * sqrt (7/pi)     * (5 * (cos theta)^3 - 3 * cos theta)
p32 (theta, phi) = (1/4)  * sqrt (105/4*pi) * (cos theta * (sin theta)^2 * cos (2*phi))

testPWithF21 :: Double -> Bool
testPWithF21 x = let
  t1 = (calcAssLegenSlow 3 0 x)    - 0.5*(5*x^3-3*x)
  t2 = (calcAssLegenSlow 4 4 x)    - 105*(1-x^2)^2
  t3 = (calcAssLegenSlow 4 (-4) x) - (105*(1-x^2)^2)/40320
  t4 = (calcAssLegenSlow 3 3 x)    - (-15*(1-x^2)**(3/2))
  t5 = (calcAssLegenSlow 3 (-3) x) - (-15*(1-x^2)**(3/2))/(-720)
  in and $ map ((< 10-8) . abs) [t1, t2, t3, t4, t5]

light :: (Double, Double) -> Double
light (theta, phi) = max 0 ((5 * cos theta) - 4) +
                     max 0 (-4 * sin (theta-pi) * cos (phi-2.5)-3)

testSH :: IO ()
testSH = let
  c   = findRealSHCoefWeight 20 (V.fromList [(10, (0, 0)), (10, (pi/4, pi/4)), (10, (pi/2, pi/2))])
  vtk = renderSphereVTK (evalRealSH c)
  in writeQuater "ODF-test-0-1" vtk

testHSH :: IO ()
testHSH = let
  c   = findRealHSHCoefWeight 20 (V.fromList [(10, (pi/2, pi/2, pi/3)), (10, (pi/4, pi/4, pi/3))])
  --vtk = renderHyperSphereVTK (evalRealHSH c)
  vtk = renderSolidSphereVTK (evalRealHSH c)
  in writeUniVTKfile ("/home/edgar/Desktop/ODF-test.vtu") False vtk

-- ====================================== Plot Space =====================================

sphere :: Int -> Int -> (Vector (Double, Double), Vector Vec3, Vector SphereCell)
sphere nPhi nTheta
  | nPhi   < 3 = sphere 3 nTheta
  | nTheta < 3 = sphere nPhi 3
  | otherwise  = (ql, V.map toCart ql, mkMesh ll)
  where
    -- phi   -> azimuthal angle
    -- theta -> polar     angle
    step_phi   = 2 * pi / (fromIntegral nPhi)
    step_theta = pi / (fromIntegral nTheta)
    funcPHI   = (step_phi *) . fromIntegral
    funcTHETA = (step_theta *) . fromIntegral
    qlm = V.fromList [ (funcTHETA nt, funcPHI np) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
    q0  = (0, 0)
    qn  = (pi, 2*pi)
    ql  = q0 `V.cons` qlm `V.snoc` qn
    toCart (theta, phi) = let
      x = sin theta * cos phi
      y = sin theta * sin phi
      z = cos theta
      in Vec3 x y z

    ll0 = V.replicate nPhi 0
    lln = V.replicate nPhi (nPhi * (nTheta-1) + 1)
    llm = V.fromList [ V.enumFromN ((n-1) * nPhi + 1) nPhi | n <- [1..nTheta-1] ]
    ll  = ll0 `V.cons` llm `V.snoc` lln

    mkMesh :: Vector (Vector Int) -> Vector SphereCell
    mkMesh vl = V.ifoldl' (\acc i l -> acc V.++ mkStrip l (vl V.! (i+1))) V.empty (V.init vl)
    mkStrip l1 l2 = let
      func i a b = SphereCell (b, l2 V.! (i+1), l1 V.! (i+1), a)
      merge      = SphereCell (V.head l2, V.head l1, V.last l1, V.last l2)
      in merge `V.cons` V.izipWith func (V.init l1) (V.init l2)

-- local instance to avoid conflict when exported.
newtype SphereCell = SphereCell (Int, Int, Int, Int)

instance RenderCell SphereCell where
  makeCell (SphereCell (a, b, c, d)) = U.fromList [a, b, c, d]
  getType _                          = VTK_QUAD

solidSphere :: Int -> Int -> Int ->
               (Vector (Double, Double, Double), Vector Vec3, Vector SolidSphereCell)
solidSphere nPhi nTheta nOmega = let
  step_omega = 2*pi/(fromIntegral nOmega)
  ws = V.enumFromN step_omega nOmega
  (polar, ps, cells) = sphere nPhi nTheta
  findR = let k = (3/4)**(1/3) in \w -> k * (w - sin w)**(1/3)
  newPS     = V.concatMap (\w -> V.map (findR w *&) ps) ws
  spherical = V.concatMap (\w -> V.map (\(p,t) -> (w, p, t)) polar) ws
  psSize = V.length ps
  func i (SphereCell (a,b,c,d)) = let
    k1 = psSize * i
    k2 = psSize * (i + 1)
    in SolidSphereCell (a+k1, b+k1, c+k1, d+k1, a+k2, b+k2, c+k2, d+k2)
  newCells  = V.concatMap (\i -> V.map (func i) cells) $ V.enumFromN 0 (nOmega-1)
  in (spherical, newPS, newCells)

newtype SolidSphereCell = SolidSphereCell (Int, Int, Int, Int, Int, Int, Int, Int)

instance RenderCell SolidSphereCell where
  makeCell (SolidSphereCell (a, b, c, d, e, f, g, h)) =
    U.fromList [a, b, c, d, e, f, g, h]
  getType _ = VTK_HEXAHEDRON

renderSphereVTK :: ((Double, Double) -> Double) -> VTK Vec3
renderSphereVTK feval = let
  (polar, ps, quads) = sphere 30 30
  func i _ = feval (polar V.!i)
  vtk  = mkUGVTK "hypersphere" (V.convert ps) quads
  attr = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderSolidSphereVTK :: ((Double, Double, Double) -> Double) -> VTK Vec3
renderSolidSphereVTK feval = let
  (polar, ps, quads) = solidSphere 30 30 30
  func i _ = feval (polar V.!i)
  vtk  = mkUGVTK "hypersphere" (V.convert ps) quads
  attr = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderHyperSphereVTK :: ((Double, Double, Double) -> Double) -> VTK Vec3
renderHyperSphereVTK feval = let
  (polar, ps, quads) = sphere 30 30
  ws = [0, 2*pi/30 .. 2*pi]
  func omega i _ = let
    (theta, phi) = polar V.!i
    in feval (omega, theta, phi)
  foo acc w = let
    attr = mkPointAttr ("Intensity - " ++ show w) (func w)
    in addDataPoints acc attr
  vtk = mkUGVTK "hypersphere" (V.convert ps) quads
  in L.foldl' foo vtk ws

renderPoints :: [Quaternion] -> VTK Vec3
renderPoints lq = let
  (pos, omega) = V.unzip . V.map (axisAngle . fromQuaternion) $ V.fromList lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  vtk  = mkUGVTK "samples" (V.convert pos) pids
  attr = mkPointAttr ("Omegas") (\a _ -> (omega) V.! a)
  in addDataPoints vtk attr

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") True
