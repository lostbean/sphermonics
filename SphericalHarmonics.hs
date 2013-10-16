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
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector         (Vector)
import           System.Random       (mkStdGen, randomRs)

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

calcK :: LMN -> Double
calcK (L l, M m, _) = let
  a = log $ fromIntegral $ 2 * l + 1
  b = logFact (l - abs m)
  c = logFact (l + abs m)
  in sqrt $ exp (a + b - log (4*pi) - c)

calcY :: LMN -> Pyramid Double -> Pyramid Double -> (Double, Double) -> Double
calcY lmn@(pl, pm, pn) pyK pyP (_, phi)
  | pm > 0    = sqrt2 * (pyK %! lmn) * cos ( m * phi) * pyP %! lmn
  | pm < 0    = sqrt2 * (pyK %! lmn) * sin (-m * phi) * pyP %! (pl, -pm, pn)
  | otherwise = let l0n = (pl, 0, pn) in (pyK %! l0n) * (pyP %!l0n)
  where m = fromIntegral (unM pm)

evalRealSH :: Pyramid Double -> (Double, Double) -> Double
evalRealSH pyC x@(theta, _) = let
  psc = pyramidStruct pyC
  l  = lmax psc
  p  = genAssLegenPyramid l (cos theta)
  k  = generatePyramid calcK psc
  ps = genPyStructLinSeq psc
  func acc lmn = acc + (pyC %! lmn) * (calcY lmn k p x)
  in L.foldl' func 0 ps

findRealSHCoef :: Int -> (Double, Double) -> Pyramid Double
findRealSHCoef l x@(theta, _) = let
  s = mkPyramidStruct (L l) FullRange NoRange
  p = genAssLegenPyramid (L l) (cos theta)
  k = generatePyramid calcK s
  func lmn = calcY lmn k p x
  in generatePyramid func s

findRealSHCoefWeight :: Int -> Vector (Double, (Double, Double)) -> Pyramid Double
findRealSHCoefWeight l xs = let
  s = mkPyramidStruct (L l) FullRange NoRange
  k = generatePyramid calcK s
  n = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo (w, x@(theta, _)) = let
    p = genAssLegenPyramid (L l) (cos theta)
    in generatePyramid (\lmn -> w * calcY lmn k p x) s
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/n) total

findRealSHCoefWith :: Int -> ((Double, Double) -> Double) -> Vector (Double, Double) -> Pyramid Double
findRealSHCoefWith l func xs = let
  s = mkPyramidStruct (L l) FullRange NoRange
  k = generatePyramid calcK s
  n = fromIntegral $ V.length xs
  x0 = V.head xs
  xt = V.tail xs
  foo x@(theta, _) = let
    p = genAssLegenPyramid (L l) (cos theta)
    f = func x
    in generatePyramid (\lmn -> f * calcY lmn k p x) s
  total = V.foldl (\acc x -> zipPyramidWith (+) acc (foo x)) (foo x0) xt
  in mapPyramid (/n) total

-- =============================== HyperSpherical Harmonics ==============================

calcHyperK :: LMN -> Double
calcHyperK (L l, M m, N n) = let
  a = log $ fromIntegral $ 2 * l + 1
  b = logFact (l - abs m)
  c = logFact (l + abs m)
  d = logFact (n - l)
  e = logFact (n + l + 1)
  r = sqrt $ exp (a + b + log (fromIntegral $ n+1) + d - c - e)
  in 2^l * (fromIntegral $ fact l) * r / pi

calcZ :: LMN -> Pyramid Double -> Pyramid Double ->
         Pyramid Double -> (Double, Double, Double) -> Double
calcZ lmn@(pl, pm, pn) pyK pyP pyC (omega, theta, phi)
  | pm > 0    = sign * common * cos (m * phi)
  | pm < 0    = sign * common * sin (m * phi)
  | otherwise = common
  where
    sign   = if even (unM pm) then 1 else -1
    common = (pyK %! lmn) * (sin (omega/2))^l * pyC %! (pl+1, 0, pn-(N l)) * pyP %! lmn
    l0n    = (pl, 0, pn)
    m      = fromIntegral (unM pm)
    l      = fromIntegral (unL pl)

-- ================================= Test functions ======================================

testPFunc :: LMN -> (Double, Double) -> Double
testPFunc lmn@(l,_,_) x@(theta, _) = let
  s = mkPyramidStruct l FullRange NoRange
  k = generatePyramid calcK s
  p = genAssLegenPyramid l (cos theta)
  in calcY lmn k p x

p70 (theta, phi) = (1/32) * sqrt (15/pi)    * (429 * (cos theta)^7 - 693 * (cos theta)^5 + 315 * (cos theta)^3 - 35 * cos theta)
p40 (theta, phi) = (3/16) * sqrt (1/pi)     * (35 * (cos theta)^4 - 30 * (cos theta)^2 + 3)
p30 (theta, phi) = (1/4)  * sqrt (7/pi)     * (5 * (cos theta)^3 - 3 * cos theta)
p32 (theta, phi) = (1/4)  * sqrt (105/4*pi) * (cos theta * (sin theta)^2 * cos (2*phi))

renderSpherePVTK :: LMN -> VTK Vec3
renderSpherePVTK lmn = let
  (polar, ps, quads) = sphere 30 30
  func i _ = testPFunc lmn (polar V.!i)
  vtk  = mkUGVTK "hypersphere" (V.convert ps) quads
  attr = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

testPWithF21 :: Double -> Bool
testPWithF21 x = let
  t1 = (calcAssLegenSlow 3 0 x)    - 0.5*(5*x^3-3*x)
  t2 = (calcAssLegenSlow 4 4 x)    - 105*(1-x^2)^2
  t3 = (calcAssLegenSlow 4 (-4) x) - (105*(1-x^2)^2)/40320
  t4 = (calcAssLegenSlow 3 3 x)    - (-15*(1-x^2)**(3/2))
  t5 = (calcAssLegenSlow 3 (-3) x) - (-15*(1-x^2)**(3/2))/(-720)
  in and $ map ((< 10-8) . abs) [t1, t2, t3, t4, t5]

testIndexAccess :: PyramidStructure -> Bool
testIndexAccess s = let
  size = linSize s
  func i = i == getPos s (dbg ">>" $ getLMN s i)
  t1 = and $ map func [0 .. size - 1]
  t2 = and $ zipWith (\i x -> getLMN s i == x) [0 .. size - 1] (genPyStructLinSeq s)
  in t1 && t2

testN :: Int -> [(Double, Double)]
testN n = let
  g1 = mkStdGen 666
  g2 = mkStdGen 667
  theta = randomRs (0, pi) g1
  phi   = randomRs (0, 2*pi) g2
  in take n $ zip theta phi

light :: (Double, Double) -> Double
light (theta, phi) = max 0 ((5 * cos theta) - 4) +
                     max 0 (-4 * sin (theta-pi) * cos (phi-2.5)-3)

testSH :: IO ()
testSH = let
  xs  = testN 10000
  c   = findRealSHCoefWeight 10 (V.fromList [(10, (0, 0)), (30, (pi/2, 0)), (10, (0, pi/2))])
  vtk = renderSphereVTK c
  in writeQuater "ODF-test-0-0" vtk

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

renderSphereVTK :: Pyramid Double -> VTK Vec3
renderSphereVTK cCoef = let
  (polar, ps, quads) = sphere 30 30
  func i _ = evalRealSH cCoef (polar V.!i)
  vtk  = mkUGVTK "hypersphere" (V.convert ps) quads
  attr = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderPoints :: [Quaternion] -> VTK Vec3
renderPoints lq = let
  (pos, omega) = V.unzip . V.map (axisAngle . fromQuaternion) $ V.fromList lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  vtk  = mkUGVTK "samples" (V.convert pos) pids
  attr = mkPointAttr ("Omegas") (\a _ -> (omega) V.! a)
  in addDataPoints vtk attr

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") True
