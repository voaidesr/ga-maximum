module Main where

import System.Random
import Text.Printf (printf)

data Config = Config
  { popSize :: Int,
    domain :: (Double, Double),
    coeffs :: (Double, Double, Double),
    precision :: Int,
    crossProb :: Double,
    mutProb :: Double,
    numSteps :: Double,
    chromLen :: Int
  }
  deriving (Show)

type Bit = Int

type Chrom = [Bit]

data Individual = Individual
  { chrom :: Chrom,
    trueVal :: Double,
    fitness :: Double
  }
  deriving (Show, Eq)

-- math helpers

chromLength :: (Double, Double) -> Int -> Int
chromLength (a, b) p = ceiling . logBase 2 $ (b - a) * 10 ^^ p

decode :: Config -> Chrom -> Double
decode = undefined

main :: IO ()
main = do
  -- Function: -x^2 + x + 2, Domain: [-1, 2], Precision: 6
  let testConfig =
        Config
          { popSize = 20,
            domain = (-1.0, 2.0),
            coeffs = (-1.0, 1.0, 2.0),
            precision = 6,
            crossProb = 0.25,
            mutProb = 0.01,
            numSteps = 50,
            chromLen = 0 -- We will calculate this next
          }

  let l = chromLength (domain testConfig) (precision testConfig)

  printf "Testing Chromosome Length:\n"
  printf "Domain: %s, Precision: %d\n" (show $ domain testConfig) (precision testConfig)
  printf "Calculated Length (L): %d bits\n\n" l
