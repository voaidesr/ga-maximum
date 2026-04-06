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
_bitsToInt :: [Bit] -> Int
_bitsToInt = foldl (\acc bit -> 2 * acc + bit) 0

chromLength :: (Double, Double) -> Int -> Int
chromLength (a, b) p = ceiling . logBase 2 $ (b - a) * 10 ^^ p

decode :: Config -> Chrom -> Double
decode cfg chr =
  let (a, b) = domain cfg
      l = chromLen cfg
      val = fromIntegral . _bitsToInt $ chr
   in a + (val * (b - a) / (2 ^ l - 1))

calcFitness :: Config -> Double -> Double
calcFitness cfg val =
  let (a, b, c) = coeffs cfg
   in a * val ** 2 + b * val + c

makeIndividual :: Config -> Chrom -> Individual
makeIndividual cfg chr =
  let val = decode cfg chr
      fit = calcFitness cfg val
   in Individual chr val fit

main :: IO ()
main = do
  -- Function: -x^2 + x + 2, Domain: [-1, 2], Precision: 6
  let d = (-2.0, 2.0)
  let p = 3
  let testConfig =
        Config
          { popSize = 20,
            domain = d,
            coeffs = (-1.0, 1.0, 2.0),
            precision = p,
            crossProb = 0.25,
            mutProb = 0.01,
            numSteps = 50,
            chromLen = chromLength d p
          }

  let l = chromLength (domain testConfig) (precision testConfig)

  printf "Testing Chromosome Length:\n"
  printf "Domain: %s, Precision: %d\n" (show $ domain testConfig) (precision testConfig)
  printf "Calculated Length (L): %d bits\n\n" l

  printf "%d\n" (_bitsToInt [1, 0, 1, 0])

  printf "%d\n" (chromLen testConfig)

  printf "%f\n" (decode testConfig [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1])
