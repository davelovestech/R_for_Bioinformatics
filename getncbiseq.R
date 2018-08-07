# I am working on the 'Retrieving genome sequence data using SeqinR' section from 
# page 15 (PDF page 19) of 'A Little Book of R for Bioinformatics, Release 0.1'.

# This is the error message I get when I run the code (below)
# Error in socketConnection(host = host, port = port, server = server, blocking = blocking,  : 
# cannot open the connection
# In addition: Warning message:
# In socketConnection(host = host, port = port, server = server, blocking = blocking,  :
# pbil.univ-lyon1.fr:5558 cannot be opened
# Error in choosebank(db, timeout = 120) : 
# I wasn't able to open the socket connection:
#   o Check that your are connected to the internet.
#   o Check that port 5558 is not closed by a firewall.
#   o Try to increase timeout value (current is 5 seconds).

# My OS is Ubuntu 16.04. I tried updating the original timeout=5 to timeout=120.
# I tried 'sudo ufw allow 5558' and no change.  
# tried 'sudo iptables -A INPUT -m state --state NEW -m tcp -p tcp --dport 5558 -j ACCEPT

getncbiseq <- function(accession)
{
  require("seqinr") # this function requires the SeqinR R package
  # first find which ACNUC database the accession is stored in:
  dbs <- c("genbank","regseq","refseqViruses","bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs)
  {
    db <- dbs[i]
    choosebank(db, timeout=120)
    # check if the sequence is in ACNUC database 'db':
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE)
    if (!(inherits(resquery, "try-error")))
    {
      queryname <- "query2"
      thequery <- paste("AC=",accession,sep="")
      query('queryname','thequery')
      # see if a sequence was retrieved:
      seq <- getSequence(query2$req[[1]])
      closebank()
      return(seq)
    }
    closebank()
  }
  print(paste("ERROR: accession",accession,"was not found"))
}
dengueseq <- getncbiseq("NC_001477")
