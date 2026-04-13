module Main where

import Data.Function (on)
import Data.List
import qualified Data.Vector as V
import System.IO
import System.Random
import Text.Printf (hPrintf)

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

selectionFit :: Config -> [Individual] -> StdGen -> ([Individual], StdGen)
selectionFit cfg pop rng =
  let probs = calcProb pop
      q = V.fromList $ calcCumulative probs
      popV = V.fromList pop
      size = popSize cfg
      (randomN, rng_) = genNRandomDoubles size rng
      selected = map (\u -> popV V.! binarySearch q u) randomN
   in (selected, rng_)

selectionElite :: [Individual] -> Individual
selectionElite = maximumBy (compare `on` fitness)

_groupParents :: [a] -> ([(a, a)], [a])
_groupParents [] = ([], [])
_groupParents [x] = ([], [x])
_groupParents (x : y : xs) =
  let (pairs, left) = _groupParents xs
   in ((x, y) : pairs, left)

_chooseCross :: Config -> [Individual] -> StdGen -> ([Individual], [Individual], StdGen)
_chooseCross cfg pop rng =
  let p = crossProb cfg
      (toss, rng_) = genNRandomDoubles (length pop) rng
      (p1_, p2_) = partition (\(_, prob) -> prob < p) $ zip pop toss
      p1 = map fst p1_
      p2 = map fst p2_
   in (p1, p2, rng_)

_crossOver :: Config -> (Individual, Individual) -> StdGen -> ([Individual], StdGen)
_crossOver cfg (i1, i2) rng =
  let c1 = chrom i1
      c2 = chrom i2
      l = chromLen cfg
      (k, rng_) = randomR (1 :: Int, l - 1) rng -- could set to 0, l if I wanted not to allow crossover to happen sometimes
      (c11, c12) = splitAt k c1
      (c21, c22) = splitAt k c2
      i1_ = makeIndividual cfg (c11 ++ c22)
      i2_ = makeIndividual cfg (c21 ++ c12)
   in ([i1_, i2_], rng_)

crossOver :: Config -> [(Individual, Individual)] -> StdGen -> ([Individual], StdGen)
crossOver _ [] rng = ([], rng)
crossOver cfg (x : xs) rng =
  let (l1, rng1) = _crossOver cfg x rng
      (l2, rng2) = crossOver cfg xs rng1
   in (l1 ++ l2, rng2)

complementBit :: Bit -> Bit
complementBit 0 = 1
complementBit 1 = 0
complementBit _ = 0

mutateC :: Config -> Chrom -> StdGen -> (Chrom, StdGen)
mutateC cfg chrm rng =
  let l = chromLen cfg
      mutP = mutProb cfg
      (probs, rng_) = genNRandomDoubles l rng
      chromProbs = zip chrm probs
      newChrom = map (\(x, p) -> if p < mutP then complementBit x else x) chromProbs
   in (newChrom, rng_)

mutateI :: Config -> Individual -> StdGen -> (Individual, StdGen)
mutateI cfg ind rng =
  let chrm = chrom ind
      (chrm_, rng_) = mutateC cfg chrm rng
      ind_ = makeIndividual cfg chrm_
   in (ind_, rng_)

mutateIs :: Config -> [Individual] -> StdGen -> ([Individual], StdGen)
mutateIs _ [] rng = ([], rng)
mutateIs cfg (x : xs) rng =
  let (x_, rng1) = mutateI cfg x rng
      (xs_, rng2) = mutateIs cfg xs rng1
   in (x_ : xs_, rng2)

formatChrom :: Chrom -> String
formatChrom = concatMap show

logMeta :: Handle -> Config -> IO ()
logMeta h cfg = do
  let (a, b, c) = coeffs cfg
  let (s, d) = domain cfg
  hPrintf
    h
    "Polynomial: %.2f X^2 + %.2f X + %.2f\n"
    a
    b
    c
  hPrintf
    h
    "Search interval [%.2f, %.2f]\n\n"
    s
    d

logI :: Handle -> Individual -> IO ()
logI h ind = do
  hPrintf
    h
    "Chrom: %s | Val: %7.4f | Fit: %7.4f\n"
    (formatChrom $ chrom ind)
    (trueVal ind)
    (fitness ind)

-- logP :: Handle -> [Individual] -> IO ()
-- logI h pop = do
-- let probs = calcProb pop

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

  rng <- getStdGen
  let (initialPop, _) = makePopulation testConfig rng

  withFile "Evolutie.txt" WriteMode $ \h -> do
    hPutStrLn h "Start Algorithm:\n"
    logMeta h testConfig
    mapM_ (logI h) initialPop
