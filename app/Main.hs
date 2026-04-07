module Main where

import qualified Data.Vector as V
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

genNRandomBits :: Int -> StdGen -> ([Bit], StdGen)
genNRandomBits 0 rng = ([], rng)
genNRandomBits n rng =
  let (bit, rng1) = randomR (0 :: Int, 1) rng
      (bits, rng2) = genNRandomBits (n - 1) rng1
   in (bit : bits, rng2)

genNRandomDoubles :: Int -> StdGen -> ([Double], StdGen)
genNRandomDoubles 0 rng = ([], rng)
genNRandomDoubles n rng =
  let (x, rng1) = randomR (0 :: Double, 1) rng
      (xs, rng2) = genNRandomDoubles (n - 1) rng1
   in (x : xs, rng2)

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

_makePopulation :: Config -> Int -> StdGen -> ([Individual], StdGen)
_makePopulation _ 0 rng = ([], rng)
_makePopulation cfg n rng =
  let (individual, rng1) = makeRandomIndividual cfg rng
      (individuals, rng2) = _makePopulation cfg (n - 1) rng1
   in (individual : individuals, rng2)

makePopulation :: Config -> StdGen -> ([Individual], StdGen)
makePopulation cfg = _makePopulation cfg (popSize cfg)

-- selection probability based on fitness
calcProb :: [Individual] -> [Double]
calcProb population =
  let totalFit = sum (map fitness population)
   in map (\x -> fitness x / totalFit) population

calcCumulative :: [Double] -> [Double]
calcCumulative = scanl (+) 0.0

binarySearch :: V.Vector Double -> Double -> Int
binarySearch vect u = search 0 (V.length vect - 2)
  where
    search low high
      | low >= high = low
      | u < vect V.! mid = search low (mid - 1)
      | u >= vect V.! (mid + 1) = search (mid + 1) high
      | otherwise = mid
      where
        mid = (low + high) `div` 2

selection :: Config -> [Individual] -> StdGen -> ([Individual], StdGen)
selection cfg pop rng =
  let probs = calcProb pop
      q = V.fromList $ calcCumulative probs
      size = popSize cfg
      (randomN, rng_) = genNRandomDoubles size rng
      selected = map (\u -> pop !! binarySearch q u) randomN
   in (selected, rng_)

main :: IO ()
main = do
  let d = (-1.0, 2.0)
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

  rng <- getStdGen
  let (initialPop, _) = makePopulation testConfig rng

  printf "Testing Chromosome Length:\n"
  printf "Domain: %s, Precision: %d\n" (show $ domain testConfig) (precision testConfig)
  printf "Calculated Length (L): %d bits\n\n" l

  printf "%d\n" (_bitsToInt [1, 0, 1, 0])

  printf "%d\n" (chromLen testConfig)

  printf "%f\n" (decode testConfig [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1])

  -- mapM_ print initialPop
  -- mapM_ print (calcProb initialPop)
  --
  -- mapM_ print (calcCumulative (calcProb initialPop))
  let (selected, _) = selection testConfig initialPop rng

  printf "SELECTED:\n"

  mapM_ print selected
