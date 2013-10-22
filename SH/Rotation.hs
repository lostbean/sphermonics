{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Hammer.Texture.SH.Rotation
       (rotM ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Packed.Matrix  as M

import           Data.Vector         (Vector)
import           Data.Packed.Matrix  ((><))

-- import           Hammer.Math.Algebra
import           Hammer.Texture.SH.Pyramid
import           Hammer.Texture.SH.SupportFunctions
import           Hammer.Texture.SphericalHarmonics

import           Debug.Trace
dbg s x = trace (show s L.++ show x) x



-- | SH Rotation matrix for polar-angle rotations. PhD Thesis J.K. Mason eq (52)
rotM :: Int -> (Double, Double, Double) -> Vector (M.Matrix Double)
rotM llen r@(omega, theta, phi) = let
  z = genZPyramid (2*llen) r

  func l (m', m) (L lb, MF mu) = let
    ld  = fromIntegral l
    lbd = fromIntegral lb
    k1  = (pi / (2 * ld + 1)) * sqrt (2 * (2 * lbd + 1))
    -- calcClebshGordan j1 m1 j2 m2 j m
    c1 = calcClebshGordan (l, m') (lb, mu) (l, m)
    z1 = z %! (N (2*l), L lb, MF mu)
    in k1 * c1 * z1

  calcAtM l (m', m) = let
    ps = genLinSeq (2*llen)
    in L.foldl' (\acc p -> acc + func l (m', m) p) 0 ps
  mL l = ((2*l+1) >< (2*l+1)) $ map (calcAtM l) [(i, j) | i <- [-l..l], j <- [-l..l]]
  in V.generate (llen+1) mL
