{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE DeriveGeneric              #-}

module Texture.SH.Pyramid
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
         , maxKey
         )
       , getLinSize
       , getMaxStack
       , generatePyramid
       , unsafeMkPyramid
       , zipPyramidWith
       , mapPyramid
       , (%!)
       ) where

import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           GHC.Generics        (Generic)

import           Data.Binary
import           Data.Vector.Binary ()

--import           Debug.Trace

class PyraKey k where
  isKeyInRange  :: k -> Bool
  getKeyPos     :: k -> Int
  getStackLevel :: k -> Int
  getMaxKey     :: Int -> k
  genLinSeq     :: Int -> [k]

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
    in (sumN l) + l - m
  getStackLevel = unL . fst
  getMaxKey l = (L l, M 0)
  genLinSeq l = [ (li, mi)
                | li <- [0 .. L l]
                , mi <- let m = M (unL li) in [m, (m-1) .. 0]
                ]

instance PyraKey (L, MF) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (L l, MF m) = l >= 0 && abs m <= l
  getKeyPos    (L l, MF m) = let
    -- The stack size along L is given by (2l+1)^2
    m' = l - m
    -- find the cumulative size of until the previous stack
    -- therefore it is sum (2(l-1)+1) = sum (2l-1)
    -- that turns to be: 2*sumN l - sum1 l
    mStack = (l+1) * (l-1) + 1
    in mStack + m'
  getStackLevel = unL . fst
  getMaxKey l = (L l, MF (-l))
  genLinSeq l = [ (li, mi)
                | li <- [L 0 .. L l]
                , mi <- let m = unL li in [MF m, MF (m-1) .. MF (-m)]
                ]

instance PyraKey (N, L) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (N n, L l) = n >= 0 && l >= 0 && l <= n
  getKeyPos    (N n, L l) = sumN n + l
  getStackLevel = unN . fst
  getMaxKey n = (N n, L n)
  genLinSeq n = [ (ni, li)
                | ni <- [N 0 .. N n]
                , li <- let l = unN ni in [L 0 .. L l]
                ]

instance PyraKey (N, L, M) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (N n, L l, M m) = n >= 0 && l >= 0 && l <= n && m >= 0 && m <= l
  getKeyPos    (N n, L l, M m) = let
    -- find the cumulative size of until the previous stack
    -- for n>= 0 then nlmStack = [1, 4, 10, 20 ..]
    -- for l>= 0 then nlStack  = [1, 3, 6, 10  ..]
    in sumHexPyramid (n `quot` 2) + sumN l + (l - m)
  getStackLevel (n, _, _) = unN n
  getMaxKey i = let n = if even i then i else i - 1 in (N n, L n, M 0)
  genLinSeq n = [ (ni, li, mi)
                | ni <- [0, 2 .. N n]
                , li <- let l = unN ni in [0 .. L l]
                , mi <- let m = unL li in [M m, M (m-1) .. M 0]
                ]

instance PyraKey (N, L, MF) where
  {-# INLINE isKeyInRange #-}
  {-# INLINE getKeyPos    #-}
  {-# INLINE getMaxKey    #-}
  {-# INLINE genLinSeq    #-}
  isKeyInRange (N n, L l, MF m) = n >= 0 && l >= 0 && l <= n && abs m <= l
  getKeyPos    (N n, L l, MF m) = let
    -- find the cumulative size of until the previous stack
    -- for n>= 0 then nlStack = [1, 9, 25 ..] ((2i+1)²)
    -- where i = n `quot` n
    -- for l>= 0 then lStack  = [1, 3, 5, 7, 9 ..] (2l+1)
    m' = l - m
    -- (2(i-1) + 1)^2 == (2i - 1)^2
    -- 4i^2 - 4i + 1
    i = n `quot` 2
    nlStack = 4 * sumN2 i - 4 * sumN i + sum1 i - 1
    -- normal index n=0,1,2..
    --nlStack = sumN2 n
    lStack  = l * l
    in nlStack + lStack + m'
  getStackLevel (n, _, _) = unN n
  getMaxKey i = let n = if even i then i else i - 1 in (N n, L n, MF (-n))
  genLinSeq n = [ (ni, li, mi)
                | ni <- [0, 2 .. N n]
                , li <- let l = unN ni in [0 .. L l]
                , mi <- let m = unL li in [MF m, MF (m-1) .. MF (-m)]
                ]

{-# INLINE sumHexPyramid #-}
sumHexPyramid :: Int -> Int
sumHexPyramid n = n * (n + 1) * (4*n - 1) `quot` 6

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

newtype L = L { unL :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)
newtype M = M { unM :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)
newtype N = N { unN :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)

newtype LF = LF { unLF :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)
newtype MF = MF { unMF :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)
newtype NF = NF { unNF :: Int } deriving (Show, Eq, Num, Enum, Ord, Generic)

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
  { maxKey   :: k
  , pyramid  :: U.Vector a
  } deriving (Show)

instance (Binary a, U.Unbox a, Binary k ,PyraKey k)=> Binary (Pyramid k a) where
  put (Pyramid k p) = do
    put k
    put p

  get = do
    k <- get
    p <- get
    return (Pyramid k p)

instance Binary L
instance Binary M
instance Binary N
instance Binary LF
instance Binary MF
instance Binary NF

unsafeMkPyramid :: (U.Unbox a, PyraKey k)=> Int -> U.Vector a -> Pyramid k a
unsafeMkPyramid s = Pyramid (getMaxKey s)

getLinSize :: (PyraKey k)=> k -> Int
getLinSize = (+1) . getKeyPos

getMaxStack :: (PyraKey k, U.Unbox a)=> Pyramid k a -> Int
getMaxStack = getStackLevel . maxKey

generatePyramid :: (U.Unbox a, PyraKey k)=> (k -> a) -> Int -> Pyramid k a
generatePyramid func smax = let
  k   = getMaxKey smax
  ks  = V.fromList $ genLinSeq smax
  vec = U.convert $ V.map func ks
  in Pyramid k vec

zipPyramidWith :: (U.Unbox a, U.Unbox b, U.Unbox c, PyraKey k, Eq k)=>
                  (a -> b -> c) -> Pyramid k a -> Pyramid k b -> Pyramid k c
zipPyramidWith func pya pyb
  | (maxKey pya) == (maxKey pyb) = pya {pyramid = v}
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

<<<<<<< HEAD
=======
{--
>>>>>>> kernel
-- ========================================= Test ========================================

testIndexAccess :: Int -> Bool
testIndexAccess n = let
  ks = V.fromList $ genLinSeq n :: V.Vector (N, L, M)
  func i p = trace (show (i, p, getKeyPos p)) $ i == getKeyPos p
  in V.and $ V.imap func ks
<<<<<<< HEAD
=======
--}
>>>>>>> kernel
