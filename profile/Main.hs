{-# LANGUAGE RecordWildCards #-}
module Main where

import Options.Applicative
import Control.Monad

import Texture.SphericalHarmonics


data Tester =
  Tester
  { run_rot_SH  :: Bool
  , run_rot_HSH :: Bool
  , run_sym_SH  :: Bool
  , run_sym_HSH :: Bool
  , run_fam_SH  :: Bool
  , run_fam_HSH :: Bool
  , run_fit_HSH :: Bool
  } deriving (Show)

tester :: Parser Tester
tester = Tester
  <$> switch
      (  long "rot-SH"
      <> help "Run test on SH rotation." )
  <*> switch
      (  long "rot-HSH"
      <> help "Run test on HSH rotation." )
  <*> switch
      (  long "symm-SH"
      <> help "Run test on SH symmetry." )
  <*> switch
      (  long "symm-HSH"
      <> help "Run test on HSH symmetry." )
  <*> switch
      (  long "base-SH"
      <> help "Plot SH base functions." )
  <*> switch
      (  long "base-HSH"
      <> help "Plot HSH base functions." )
  <*> switch
      (  long "fit-HSH"
      <> help "fit HSH on 1000 points" )

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> tester)
      ( fullDesc
      <> progDesc "Test and profile sledge library."
      <> header "Hammer library" )

run :: Tester -> IO ()
run Tester{..} = do
  when run_rot_SH  testRotSH
  when run_rot_HSH testRotHSH
  --when run_sym_SH  testSymmSH  -- no implemented yet
  when run_sym_HSH testSymmHSH
  when run_fam_SH  (plotSHFuncFamily 10)
  when run_fam_HSH (plotHSHFuncFamily 10)
