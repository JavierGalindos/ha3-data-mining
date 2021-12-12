# Exercise 1: Text Data Mining
# Javier Galindos

rm(list = ls()) # Remove all the objects we created so far.

library("rstudioapi")  

# Set working directory
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()  

cname <- file.path("./Texts")
dir(cname) #check loaded texts

# load package tm  which is framework for text data mining / NLP
library(tm)
library(mixtools) # EM
library(fpc)
library(cluster)
library(ggplot2)

# load available documents
docs <- VCorpus(DirSource(cname))   #corpus or set of texts
#print(summary(docs))

inspect(docs[2]) #some data about a loaded text - order 2
writeLines(as.character(docs[2])) #content of text 2

# Preprocessing
docs <- tm_map(docs,removePunctuation)   #remove punctuation symbols
docs <- tm_map(docs, removeNumbers)   #remove numbers
docs <- tm_map(docs, tolower)   #remove capitalization - to lowcase
docs <- tm_map(docs, removeWords, c("’ll", "’re","’ve"))   #remove particular words
docs <- tm_map(docs, removeWords, c("the", "and", stopwords("english")))   #remove common words such as "a, and, also, the"
docs <- tm_map(docs, stripWhitespace)  #remove whitespaces from previous eliminations
docs <- tm_map(docs, PlainTextDocument)   #prepare de document as text



#Work with the data
#we use a Document-Term Matrix (DTM) representation: documents as the rows, terms/words as the columns, frequency of the term in the document as the entries. Because the number of unique words in the corpus the dimension can be large.
dtm <- DocumentTermMatrix(docs)   #create document-term matrix
dim(dtm)
print(dtm)

tdm <- TermDocumentMatrix(docs)   #creating the transpose of dtm
print(tdm)

# Represent sentences numerically (Weighted TF-IDF)
tdm.tfidf <- tm::weightTfIdf(tdm)
print(tdm.tfidf)


#exploring data
freq <- colSums(as.matrix(dtm))   #organize terms by frequency
length(freq)

head(table(freq), 10) #check the less occuring words - top number is frequency of appearing, bottom number the number of words
tail(table(freq), 10) #check the most occurring words

dtms <- removeSparseTerms(dtm, 0.1) #remove sparse - less frequent items, matrix that is only 10% empty space
print(dtms)

freq <- colSums(as.matrix(dtms))   
freq   

# table after removing sparse terms
freq <- colSums(as.matrix(dtms))  
print(freq)

print(findFreqTerms(dtm, lowfreq=30)) #show words that appear 30 or more times

#plotting

# word frequencies

wf <- data.frame(word=names(freq), freq=freq)
p <- ggplot(subset(wf, freq>5), aes(word, freq))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
print(p)

#or ordered...
p <- ggplot(subset(wf, freq>5), aes(x = reorder(word, -freq), y = freq)) +
  geom_bar(stat = "identity", color="indianred", fill = "indianred") + 
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  ggtitle("Most common words Biden")
print(p)  

#find correlations between terms
findAssocs(dtm, c("will" , "american"), corlimit=0.85) #find correlations between words

#cloud of words
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

set.seed(1234)
freq = data.frame(sort(colSums(as.matrix(dtm)), decreasing=TRUE))
wordcloud::wordcloud(rownames(freq), freq[,1], scale=c(4, .1), max.words=50, colors=brewer.pal(6, "Dark2")) #most frequently used 50 words

#plot words that occur at least 50 times
wordcloud::wordcloud(rownames(freq), freq[,1], scale=c(4, .1), min.freq=25, colors=brewer.pal(6, "Dark2"))  #words used at least 25 times or more

  

# K-means
dtms <- removeSparseTerms(dtm, 0.15)
# d <- dist(t(dtms), method="manhattan")
#  Cosine similarity
d <- proxy::dist(as.matrix(t(dtms)) , method = "cosine")
clustering.kmeans <- kmeans(d, 3)
clusplot(as.matrix(d), clustering.kmeans$cluster, color=T, shade=T, labels=2, lines=0, main = "K-means")

# Hierarchical clustering
clustering.hierarchical <- hclust(d, method = "ward.D2")
plot(clustering.hierarchical, hang=-1) #plot

plot.new() #find number of clusters
plot(clustering.hierarchical, hang=-1)
groups <- cutree(clustering.hierarchical, k=3)   # "k=" defines the number of clusters you are using   
rect.hclust(clustering.hierarchical, k=3, border="red") # draw dendogram with red borders around the 6 clusters   


# HBDSCAN
clustering.dbscan <- dbscan::hdbscan(d, minPts = 5)
clusplot(as.matrix(d), clustering.dbscan$cluster, color=T, shade=T, labels=2, lines=0)

# PAM
clustering.pam <- pam(d, 3)
clusplot(as.matrix(d), clustering.pam$cluster, color=T, shade=T, labels=2, lines=0, main = "K-medoids")

