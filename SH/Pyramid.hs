{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Hammer.Texture.SH.Pyramid
       ( L (..)
       , M (..)
       , N (..)
       , LMN
       , Pyramid
         ( pyramid
         , pyramidStruct
         )
       , PyramidRange (..)
       , PyramidStructure
         ( lmax
         , linSize
         , mRange
         , nRange
         , mSize
         , nSize
         , layerLast
         , layerSize
         )
       , mkPyramidStruct
       , unsafeMkPyramid
       , generatePyramid
       , getPyramidLen
       , getLayerSize
       , getLMN
       , getPos
       , getMRange
       , getNRange
       , genPyStructLinSeq
       , zipPyramidWith
       , mapPyramid
       , isSamePyramidStruct
       , (%!)
       ) where

import qualified Data.Vector.Unboxed as U

-- ================================= Pyramid =====================================

newtype L = L { unL :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype M = M { unM :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype N = N { unN :: Int } deriving (Show, Eq, Num, Enum, Ord)

type LMN = (L, M, N)

data PyramidRange = FullRange
                  | HalfRange
                  | NoRange
                    deriving (Show, Eq)

data (U.Unbox a)=> Pyramid a =
  Pyramid
  { pyramidStruct :: PyramidStructure
  , pyramid       :: U.Vector a
  } deriving (Show)

data PyramidStructure =
  PyramidStructure
  { lmax      :: L
  , linSize   :: Int
  , mRange    :: PyramidRange
  , nRange    :: PyramidRange
  , mSize     :: U.Vector Int
  , nSize     :: U.Vector Int
  , layerLast :: U.Vector Int
  , layerSize :: U.Vector Int
  } deriving (Show)

mkPyramidStruct :: L -> PyramidRange -> PyramidRange -> PyramidStructure
mkPyramidStruct l@(L lm) mr nr = let
  lls = U.postscanl (+) 0 lss
  lss = U.zipWith   (*) mss nss
  mss = U.generate  (lm+1) (getPyramidLen mr . L)
  nss = U.generate  (lm+1) (getPyramidLen nr . L)
  in  PyramidStructure { lmax      = l
                       , linSize   = if lm >= 0 then U.last lls else 0
                       , mRange    = mr
                       , nRange    = nr
                       , mSize     = mss
                       , nSize     = nss
                       , layerLast = lls
                       , layerSize = lss
                       }

unsafeMkPyramid :: (U.Unbox a)=> PyramidStructure -> U.Vector a -> Pyramid a
unsafeMkPyramid = Pyramid

getPyramidLen :: PyramidRange -> (L -> Int)
getPyramidLen pyrange = case pyrange of
  FullRange -> (+1) . (* 2) . unL
  HalfRange -> (+1) . unL
  _         -> const 1

getLayerSize :: PyramidRange -> PyramidRange -> (L -> Int)
getLayerSize m n = let
  fm = getPyramidLen m
  fn = getPyramidLen n
  in \i -> fm i * fn i

getLMN :: PyramidStructure -> Int -> LMN
getLMN PyramidStructure{..} i = let
  offSet pr = case pr of
    FullRange -> \x -> x-l
    _         -> id
  l  = maybe (unL lmax) id (U.findIndex (> i) layerLast)
  ll = layerLast U.! l
  ns = nSize     U.! l
  ms = mSize     U.! l
  t  = i - (ll - ms*ns)
  (m,n) = t `quotRem` ns
  m' = offSet mRange m
  n' = offSet nRange n
  in (L l, M m', N n')

getPos :: PyramidStructure -> LMN -> Int
getPos PyramidStructure{..} (L l, M m, N n) = let
  offSet pr = case pr of
    FullRange -> \x -> x+l
    _         -> id
  ns = nSize     U.! l
  ll = layerLast U.! l
  ls = layerSize U.! l
  m' = offSet mRange m
  n' = offSet nRange n
  in (ll - ls) + m'*ns + n'

getMRange :: PyramidStructure -> L -> (M, M)
getMRange PyramidStructure{..} (L l) = let
  ms = (mSize U.! l) - 1
  in case mRange of
    FullRange -> (M $ l-ms, M $ ms-l)
    _         -> (M 0, M ms)

getNRange :: PyramidStructure -> L -> (N, N)
getNRange PyramidStructure{..} (L l) = let
  ns = (nSize U.! l) - 1
  in case nRange of
    FullRange -> (N $ l-ns, N $ ns-l)
    _         -> (N 0, N ns)

genPyStructLinSeq :: PyramidStructure -> [LMN]
genPyStructLinSeq s@PyramidStructure{..} = let
  in [(li, mi, ni)
     | li <- [0 .. lmax]
     , mi <- let (mb, mu) = getMRange s li in [mb..mu]
     , ni <- let (nb, nu) = getNRange s li in [nb..nu]]

generatePyramid :: (U.Unbox a)=> (LMN -> a) -> PyramidStructure -> Pyramid a
generatePyramid func p = let
  vec = U.generate (linSize p) (func . getLMN p)
  in Pyramid p vec

zipPyramidWith :: (U.Unbox a, U.Unbox b, U.Unbox c)=> (a -> b -> c) -> Pyramid a -> Pyramid b -> Pyramid c
zipPyramidWith func pya pyb
  | isSamePyramidStruct (pyramidStruct pya) (pyramidStruct pyb) = pya {pyramid = v}
  | otherwise = error "Can't zip two Pyramids with different structures."
  where
    v = U.zipWith func (pyramid pya) (pyramid pyb)

mapPyramid :: (U.Unbox a, U.Unbox b)=> (a -> b) -> Pyramid a -> Pyramid b
mapPyramid func py = py { pyramid = U.map func (pyramid py)}

isSamePyramidStruct :: PyramidStructure -> PyramidStructure -> Bool
isSamePyramidStruct psa psb = lmax   psa == lmax   psb &&
                              mRange psa == mRange psb &&
                              nRange psa == nRange psb

(%!) :: (U.Unbox a)=> Pyramid a -> LMN -> a
Pyramid{..} %! lmn = pyramid U.! (getPos pyramidStruct lmn)
