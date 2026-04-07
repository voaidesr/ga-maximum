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

genNRandomBits :: Bit -> StdGen -> ([Bit], StdGen)
genNRandomBits 0 rng = ([], rng)
genNRandomBits n rng =
  let (bit, rng1) = randomR (0 :: Int, 1) rng
      (bits, rng2) = genNRandomBits (n - 1) rng1
   in (bit : bits, rng2)

makeRandomChromo :: Config -> StdGen -> (Chrom, StdGen)
makeRandomChromo cfg rng =
  let l = chromLen cfg
      (chrm, rng_) = genNRandomBits l rng
   in (chrm, rng_)

makeIndividual :: Config -> Chrom -> Individual
makeIndividual cfg chr =
  let val = decode cfg chr
      fit = calcFitness cfg val
   in Individual chr val fit

makeRandomIndividual :: Config -> StdGen -> (Individual, StdGen)
makeRandomIndividual cfg rng =
  let (chrm, rng_) = makeRandomChromo cfg rng
      individual = makeIndividual cfg chrm
   in (individual, rng_)

makePopulation :: Config -> Int -> StdGen -> ([Individual], StdGen)
makePopulation _ 0 rng = ([], rng)
makePopulation cfg n rng =
  let (individual, rng1) = makeRandomIndividual cfg rng
      (individuals, rng2) = makePopulation cfg (n - 1) rng1
   in (individual : individuals, rng2)

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
