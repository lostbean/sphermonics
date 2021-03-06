{-# LANGUAGE RecordWildCards #-}
module Main where

import qualified Data.Vector as V
import qualified Data.HashMap.Strict as HM

import Options.Applicative
import Data.Binary
import Data.Complex

import Texture.TesseractGrid
import Texture.Orientation
import Texture.SphericalHarmonics
import Texture.HyperSphere


data TableGen =
  TableGen
  { gen_legendre   :: Maybe Int
  , gen_hsh_coeff  :: Maybe Int
  , outdir         :: String
  } deriving (Show)

genInfo :: Parser TableGen
genInfo = TableGen
  <$> (optional (option
      ( long     "legendre"
      <> metavar "L-level"
      <> help    "Generate series of legendre functions." )))
  <*> (optional (option
      ( long     "hsh-coeff"
      <> metavar "N-level"
      <> help    "Generate series of Spherical Harmonics coefficients." )))
  <*> (strOption
      ( short     'd'
      <> metavar "DIRPATH"
      <> help    "Output directory to save the tables." ))

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> genInfo)
      ( fullDesc
      <> progDesc "Test and profile sledge library."
      <> header "Hammer library" )

run :: TableGen -> IO ()
run TableGen{..} = do
  maybe (return ()) (genLegendre  outdir) gen_legendre
  maybe (return ()) (genHSHCoeff  outdir) gen_hsh_coeff

genHSHCoeff :: FilePath -> Int -> IO ()
genHSHCoeff dir li = let
  tess = genTesseractGrid2 10 sh
  sh :: Quaternion -> Pyramid (N, L, MF) (Complex Double)
  sh = genSHFunc li . quaternionToSO3
  in do
    putStrLn "Generating HSH table..."
    encodeFile (dir ++ "/" ++ "HSH-L=" ++ show li ++ ".data") tess
    putStrLn "Done!"

genLegendre :: FilePath -> Int -> IO ()
genLegendre dir li = let
  n  = 360
  ws = V.enumFromStepN 0 (2*pi/(fromIntegral n)) n
  ls = V.map (genAssLegenFullPyramid (L li) . cos) ws
  gs = V.map (genGegenbauerPyramid   (N li) . cos) ws
  in do
    putStrLn "Generating Legendre and Gegenbauer table..."
    encodeFile (dir ++ "/" ++ "legendre_coef.data") ls
    encodeFile (dir ++ "/" ++ "gegenbauer_coef.data") gs
    encodeFile (dir ++ "/" ++ "clebsh-n=" ++ show n ++ ".data") $ HM.toList memoClebsh
    putStrLn "Done!"

instance (Binary a, RealFloat a)=> Binary (Complex a) where
  put c = put (realPart c) >> put (imagPart c)
  get = do
    r <- get
    i <- get
    return $ r :+ i
