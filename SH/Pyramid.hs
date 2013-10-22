{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE FlexibleInstances          #-}

module Hammer.Texture.SH.Pyramid
       ( L (..)
       , M (..)
       , N (..)
       , LF (..)
       , MF (..)
       , NF (..)
       , l2lf, m2mf, n2nf
       , lf2l, mf2m, nf2n
       , PyraKey (..)
       , Pyramid
         ( pyramid
         , maxStack
         )
       , getLinSize
       , generatePyramid
       , unsafeMkPyramid
       , zipPyramidWith
       , mapPyramid
       , (%!)
       ) where

import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

class PyraKey k where
  isKeyInRange :: k -> Bool
  getKeyPos    :: k -> Int
  getMaxKey    :: Int -> k
  genLinSeq    :: Int -> [k]

instance PyraKey (L, M) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (L l, M m) = l >= 0 && m >= 0 && m <= l
  getKeyPos    (L l, M m) = let
    -- The stack size along L is given by (l+1)
    -- find the cumulative size of until the previous stack
    -- therefore it is sum (l)
    in sumN l + m
  getMaxKey l = (L l, M l)
  genLinSeq l = [ (li, mi)
                | li <- [L 0 .. L l]
                , mi <- [M 0 .. M (unL li)]]

instance PyraKey (L, MF) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (L l, MF m) = l >= 0 && abs m <= l
  getKeyPos    (L l, MF m) = let
    -- The stack size along L is given by (2l+1)^2
    m' = l + m
    -- find the cumulative size of until the previous stack
    -- therefore it is sum (2(l-1)+1) = sum (2l-1)
    -- that turns to be: 2*sumN l - sum1 l
    mStack = (l+1) * (l-1) + 1
    in mStack + m'
  getMaxKey l = (L l, MF l)
  genLinSeq l = [ (li, mi)
                | li <- [L 0 .. L l]
                , mi <- [MF (-(unL li)) .. MF (unL li)]]

instance PyraKey (N, L) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (N n, L l) = n >= 0 && l >= 0 && l <= n
  getKeyPos    (N n, L l) = sumN n + l
  getMaxKey n = (N n, L n)
  genLinSeq n = [ (ni, li)
                | ni <- [N 0 .. N n]
                , li <- [L 0 .. L (unN ni)]]

instance PyraKey (L, MF, NF) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (L l, MF m, NF n) = l >= 0 && abs m <= l && abs n <= l
  getKeyPos    (L l, MF m, NF n) = let
    -- The stack size along L is given by (2l+1)^2
    maxLen = 2 * l + 1
    l1 = l + 1
    m' = l + m
    n' = l + n
    -- find the cumulative size of until the previous stack
    -- therefore it is sum (2(l-1)+1)^2 = sum (2l-1)^2
    -- that turns to be: 4*sumN2 l - 4*sumN l + sum1 l
    mnStack = 4 * l * l1 * maxLen `quot` 6 - 2 * l * l1 + l
    mStack  = m' * maxLen
    in mnStack + mStack + n'
  getMaxKey l = (L l, MF l, NF l)
  genLinSeq l = [ (li, mi, ni)
                | li <- [L 0 .. L l]
                , mi <- [MF (-(unL li)) .. MF (unL li)]
                , ni <- [NF (-(unL li)) .. NF (unL li)]]

instance PyraKey (N, L, MF) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (N n, L l, MF m) = n >= 0 && l >= 0 && l <= n && abs m <= l
  getKeyPos    (N n, L l, MF m) = let
    -- The stack size along L is given by (n+1)²
    m' = l + m
    -- find the cumulative size of until the previous stack
    -- for n>= 0 then nlStack = [1, 4, 9 ..] ((n+1)²)
    -- for l>= 0 then lStack  = [1, 3, 5 ..] (2l+1)
    nlStack = sumN2 n
    lStack  = l * l
    in nlStack + lStack + m'
  getMaxKey n = (N n, L n, MF n)
  genLinSeq n = [ (ni, li, mi)
                | ni <- [0 .. N n]
                , li <- [0 .. L (unN ni)]
                , mi <- [MF (-(unL li)) .. MF (unL li)]]

{-# INLINE sumN2 #-}
sumN2 :: Int -> Int
sumN2 m = m * (m+1) * (2*m+1) `quot` 6

{-# INLINE sumN #-}
sumN :: Int -> Int
sumN m = m * (m+1) `quot` 2

{-# INLINE sum1 #-}
sum1 :: Int -> Int
sum1 m = m + 1

-- ================================= Pyramid =====================================

newtype L = L { unL :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype M = M { unM :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype N = N { unN :: Int } deriving (Show, Eq, Num, Enum, Ord)

newtype LF = LF { unLF :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype MF = MF { unMF :: Int } deriving (Show, Eq, Num, Enum, Ord)
newtype NF = NF { unNF :: Int } deriving (Show, Eq, Num, Enum, Ord)

lf2l :: LF -> L
lf2l (LF l) = L $ abs l

l2lf :: L -> LF
l2lf (L l) = LF l

mf2m :: MF -> M
mf2m (MF m) = M $ abs m

m2mf :: M -> MF
m2mf (M m) = MF m

nf2n :: NF -> N
nf2n (NF n) = N $ abs n

n2nf :: N -> NF
n2nf (N n) = NF n

data (U.Unbox a, PyraKey k)=> Pyramid k a =
  Pyramid
  { maxStack :: Int
  , pyramid  :: U.Vector a
  } deriving (Show)

unsafeMkPyramid :: (U.Unbox a, PyraKey k)=> Int -> U.Vector a -> Pyramid k a
unsafeMkPyramid = Pyramid

getLinSize :: (PyraKey k)=> k -> Int
getLinSize = (+1) . getKeyPos

generatePyramid :: (U.Unbox a, PyraKey k)=> (k -> a) -> Int -> Pyramid k a
generatePyramid func smax = let
  ks  = V.fromList $ genLinSeq smax
  vec = U.convert $ V.map func ks
  in Pyramid smax vec

zipPyramidWith :: (U.Unbox a, U.Unbox b, U.Unbox c, PyraKey k, Eq k)=>
                  (a -> b -> c) -> Pyramid k a -> Pyramid k b -> Pyramid k c
zipPyramidWith func pya pyb
  | (maxStack pya) == (maxStack pyb) = pya {pyramid = v}
  | otherwise = error "Can't zip two Pyramids with different structures."
  where
    v = U.zipWith func (pyramid pya) (pyramid pyb)

mapPyramid :: (U.Unbox a, U.Unbox b, PyraKey k)=> (a -> b) -> Pyramid k a -> Pyramid k b
mapPyramid func py = py { pyramid = U.map func (pyramid py) }

{-# INLINE (%!) #-}
{-# SPECIALISE INLINE (%!) :: Pyramid (N, L, MF) Double -> (N, L, MF) -> Double #-}
{-# SPECIALISE INLINE (%!) :: Pyramid (N, L) Double -> (N, L) -> Double #-}
{-# SPECIALISE INLINE (%!) :: Pyramid (L, M) Double -> (L, M) -> Double #-}
{-# SPECIALISE INLINE (%!) :: Pyramid (L, MF) Double -> (L, MF) -> Double #-}
(%!) :: (U.Unbox a, PyraKey k)=> Pyramid k a -> k -> a
Pyramid{..} %! key = pyramid U.! (getKeyPos key)



testIndexAccess :: Int -> Bool
testIndexAccess n = let
  ks = V.fromList $ genLinSeq n :: V.Vector (N, L, MF)
  func i p = i == getKeyPos p
  in V.and $ V.imap func ks
