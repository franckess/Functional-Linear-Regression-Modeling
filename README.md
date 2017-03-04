# Record Linkage in Consumer Products Data using Approximate String Matching and Clustering Methods

## Preface

This academic paper is submitted for the degree of Master of Science at the University of Minnesota - Twin Cities. The study was conducted under the supervision of Professor Xiaoou Li in the Department of Statistics, University of Minnesota, during Fall of 2016.

## Abstract

Consumer product data is an absolute necessity for any business to do any sort of analytics. This kind of data generally pertains to sales data of products, and is typically deaggregated by product and date-time. Using consumer products data, analyses can be done in numerous ways to help inform business decisions. However, before any data analysis can be conducted, preprocessing of data (and lots of it) is necessary. One of these tasks is the joining or aggregating of data from different sources. A large issue with this is the occurrence of duplicate or redundant records. The task of identifying duplicate or redundant records is known as record linkage. Record linkage is a data cleaning task of identifying records that belong to the same entity in a single data set or across multiple data sets. In this paper we focus on the application of record linkage in consumer products data.

The goal of this paper is to discuss methods to automate the record linkage process (identify the product line of a SKU). Since the names of these SKU will contain information about the product line (PL) but will never have exact duplicates, we will use the probablistic record linkage approach. In a probablistic record linkage, an estimate of the distance between different records are obtained, and groups are computed using distances. The process is divided into two steps. The first step is computing an estimate of distance between all pairs in a list of records. The method we will use is known as approximate string matching. Approximate String Matching (also known as Fuzzy String Matching) is a pattern matching algorithm that computes the degree of similartity between two strings, and produces a quantitative metric of distance that can be used to classify the strings as a match or not a match. Using approximate string matching, we will compute a distance between all pairs of records to obtain a distance (or dissimilarity) matrix. The second step is grouping records (SKU) into clusters (PL) by applying clustering methods. In clustering (or cluster analysis), observations are grouped into clusters in such a way that objects in the same cluster are more similar to each other than those in other clusters. Clustering requires an input of a distance matrix, which is computed in the first step, and outputs the grouping of clusters. To compare the results and performance of our approaches, several evaluation methods such as accuracy, purity, and F-measure will be used.

In this paper, several existing and modified methods were applied to both real data and simulated data.
