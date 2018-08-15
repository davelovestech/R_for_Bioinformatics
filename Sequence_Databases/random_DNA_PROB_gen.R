random_DNA_PROB_gen <- function(n, pA, pC, pG, pT) 
{
  # check to make sure that the sum of probabilities == 1
  # stopifnot(sum(pA + pC + pG + pT)==1)
  # create an empty character vector for the nucleotide string
  nucleotide_vector = character()
  # run the for loop from 1 until n
  for (i in 1:n)
  {
    # print("random number: ")
    # generate a random number between 0 and 1
    random_number <- runif(1, 0, 1)
    # print(random_number)
    # print("random nucleotide: ")
    # probability tree of increasing sums of the nucleotide frequency
    # this works because it starts at the lowest numbers and goes up 
    # it wouldn't work in the opposite direction (low to high) cause everything 
    # would sort into the lowest
    if (random_number <= pA) 
    {
      random_nucleotide <- "A"
    }
    else if (random_number <= (pA + pC)) 
    {
      random_nucleotide <- "C"
    }
    else if (random_number <= (pA + pC + pG)) 
    {
      random_nucleotide <- "G"
    }
    # the probability sum isn't necessary here because higher is the only option
    # this function doesn't take into account errors in user input ... future might
    # wanna change to an else going to error and having G be assigned in an else if
    else 
    {
      random_nucleotide <- "T"
    }
    # print(random_nucleotide)
    nucleotide_vector[i] = random_nucleotide
  }
  print(nucleotide_vector)
}
cat("How many nucleotides long? ")
n <- as.integer(readLines(file("stdin"), 1))
cat("What is frequency of A? ")
pA <- as.integer(readLines(file("stdin"), 1))
cat("What is frequency of C? ")
pC <- as.integer(readLines(file("stdin"), 1))
cat("What is frequency of G? ")
pG <- as.integer(readLines(file("stdin"), 1))
cat("What is frequency of T? ")
pT <- as.integer(readLines(file("stdin"), 1))
random_DNA_PROB_gen(n, pA, pC, pG, pT)
# random_DNA_PROB_gen(20, 0.28, 0.21, 0.22, 0.29)