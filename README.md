README
================

# Ecotron Lake Experiment

## Hello Fishes

Here is a little introduction !  
You will find in this documents every information you need on those data
and the script and analyses going with it.  
I suggest you read everything to understand what we are talking about
and to be aware of the choices and hypothesis I made.  
But as I know it’s often easier to read and discuss about something we
can see and play with, you can also begin by reading about the
experiment, the list of files and jump directly to the hand on section
to get directly into this data.

## Experiment

Information on the experimantal lake : [Aquacosm.pdf](Aquacosm.pdf)

## List and quick explanation of the different files

### Dataset in /data

#### Data bases listing the Tag_ids and the related information

See section [Tag data bases]() for more precision, else a quick summary
:

-   **BDD_original complet.xslx**: row dataset input from the
    experiment  
    see explanation : [BDD_original_complet.xlsx explanation]()
-   **data_final.xlsx**, cleaned version of BDD_original with
    unification of notation and elimination of incoherent data probably
    coming from experipmental input error (such as tag_id found in wrong
    lake)  
    see explanation : [data_final.xlsx explanation]()
-   **data_final_ws.xlsx** : data_final with additional weight and size
    for juveniles and 2022  
    see hypothesis and precision : [data_final_ws.xlsx explanation]()
-   **data_final_ws_norm.xlsx** : data_final with additional weight and
    size for juveniles and 2022 by using a normal distribution  
    this is the data used in my following analyses and the creation of
    usable subset (see section [Hands on the data, getting started]())  
    see hypothesis and precision on this data : [data_final_ws_norm.xlsx
    explanation]()

#### Other data

-   **Lake_treatment.xlsx** : a table containing the experimental
    characteristics of each lake : Nutrients enrichment, presence of
    perch (predator species), and a Treatment column summarizing the two
    see precision : [Lake_tretment.xlsx explanation]()
-   **capture_history.inp** : a file containing usable data for capture
    recapture analyses with mark, i.e. capture history and corresponding
    frequencies see construction and esplanation : [capture_history.inp
    explanation]()

### Scripts

-   **construction_BDD.Rmd**
-   **importation_data.R**

## Diving deeper in the data

### Data bases

These data bases present the result of the fisheries that happened in
the lake during the 5 years experiments.  
You will find a line for each fish observed, with the different
information collected (Tag_id,date, size, weight, Lake_capture,
Lake_released), and information about the experimental characteristics
of this fishery (Method, Passage, Session).

  

#### BDD_original_complet.xlsx explanation

Simply the first data base I was given.  
Not necessary for a simple data usage. Need data cleaning before
analyses.

  

#### data_final.xlsx explanation

Cleaned version of BDD_original_complet.xlsx, see the
[construction_BDD.Rmd](construction_BDD.Rmd) file to precisely see the
modification that where made. (Not necessary for a simple usage of the
data)  

##### Columns description

| Columns name        | Description                                                                                                                                                   |
|---------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Index**           | Number of this observation                                                                                                                                    |
| **Date**            | Date of observation (-dd/mm/yyyy-)                                                                                                                            |
| **Obs_status**      | Action done with this fish -introduction- -capture- -recapture-                                                                                               |
| **Lake_capture**    | The lake in which the fish was extracted, -NA- if it is a new introduction                                                                                    |
| **Method_capture**  | How was the fish extracted ? -Trawl- (in november each year), -Hoop net- (for the perch during summer), or -Draining- (for 2022, when the lake where emptied) |
| **Session_capture** | Only one -NA- exept in 2019 as a problem of efficiency made 2 fishry necessary (one normal -A- and a second with brushing -B-)                                |
| **Passage_net**     | Every fishery had 3 net passage -1- -2- and -3-                                                                                                               |
| **Tag_id**          | Tag identifying individually each fish bigger than 8g                                                                                                         |
| **Species**         | The species of the fish (-Pike-, -Gardon-, -Perch-, -Able-, -Goujon-)                                                                                         |
| **Weight**          | Weight (in g)                                                                                                                                                 |
| **Size**            | Size of the fish from head to the fork of the tail (in mm)                                                                                                    |
| **Lake_released**   | The lake in which the fish was released, -NA- if it is a new introduction                                                                                     |
| **Comment**         | Any additional information                                                                                                                                    |
|                     |                                                                                                                                                               |

Let’s see how it looks  
  

``` r
library(readxl)
BDD_f <- read_excel("data/data_final.xlsx", 
                    col_types = c("numeric","date","text","text","text","text","text",
                                  "text","text","numeric","numeric","text","text"), na = "")
knitr::kable(head(BDD_f))
```

| Index | Date       | Lake_capture | Methode_capture | Session_capture | Passage_net | Obs_status   | Tag_id    | Species | Weight | Size | Lake_released | Comment_obs |
|------:|:-----------|:-------------|:----------------|:----------------|:------------|:-------------|:----------|:--------|-------:|-----:|:--------------|:------------|
|     1 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274142 | gardon  |     31 |  129 | 16            | NA          |
|     2 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274145 | gardon  |     11 |   97 | 16            | NA          |
|     3 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274148 | gardon  |     11 |   98 | 1             | NA          |
|     4 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274149 | gardon  |     11 |  100 | 14            | NA          |
|     5 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274150 | gardon  |     13 |   98 | 1             | NA          |
|     6 | 2016-12-06 | NA           | NA              | NA              | NA          | introduction | 403274151 | gardon  |     10 |   88 | 3             | NA          |

  

#### data_final_ws.xlsx explanation

Same as above but weight and size where added as such

**For the juveniles** :  
Adding the mean weight and size of the measured and weighted juveniles
from the same Lake and Year

**For jan 2022** :  
The hypothesis is that fishes weight and size didn’t change between nov
2021 and jan 2022, because of the time gap (4 months only) and the fact
that winter is a slow growth period  
This can be verified because some fishes where weighted at both
occasion. This difference didn’t seem significant, I then accepted the
first hypothesis  
Weight and size where added in 2022 by taking the associated weight and
size from nov 2021

  

#### data_final_ws_norm.xlsx explanation

Exactly the same as above but

**For juveniles** : weight and size where added using a normal
distribution coming from the measured and weighted juveniles from the
same Lake and Year  
This was done for two reasons :

-   It’s closer to the biologic reality
-   It diminishes the bias toward one over-represented value

### Other useful data

#### Lake_treatment.xlsx explanation

A short table giving the characteristics of each Lake :

-   **Nutrients** : TRUE if enriched
-   **Perch** : TRUE if there is no fisheries pressure on perches,
    i.e. presence of predators
-   **Treatment** :
    -   1 for Nutrients and no Perch
    -   2 for Nutrients and Perch
    -   3 for no Nutrients and no Perch
    -   4 for no Nutrients and Perch

``` r
Lake_treatment <- read_excel("data/Lake_treatment.xlsx", col_types =c("text","logical","logical","text"))
head(Lake_treatment)
#> # A tibble: 6 x 4
#>   Lake  Nutrients Perch Treatment
#>   <chr> <lgl>     <lgl> <chr>    
#> 1 1     TRUE      FALSE 1        
#> 2 2     FALSE     FALSE 3        
#> 3 3     TRUE      TRUE  2        
#> 4 4     FALSE     TRUE  4        
#> 5 5     FALSE     TRUE  4        
#> 6 6     TRUE      TRUE  2
```

  

#### capture_history.inp explanation

This data is useful for capture recapture study and can be used with the
RMark library  
Can be used for exemple with the Jolly-Seber models from oliviergimenez
available on [one of his github
repositories](https://github.com/oliviergimenez/bayes-multistate-jollyseber)  
  
**How was it made ?**  
Construction with capture_history.R ([see the code for
details](capture_history.R)) from the data of
[data_ws_norm.xlsx](data/data_ws_norm.xlsx) after importation (thanks to
[importation_data.R](data/importation_data.R))  
Utilization of a subset of the data containing only tagged fishes from
the annual fisheries of November. Then summarizing the data to get :  

``` r
library(RMark)
#popan <- convert.inp("data/capture_history.inp")
#colnames(popan) <- c("capture_history", "frequency") # simply for better comprehension here
#head(popan)
```

## Hands on the data, getting started, quick guide

This is a little explanation of how to get started with this data  
It’s often complicated to dive into such analyses so I wanted to give
anyone who wants to use in a way in.  
There are probably other ways to look at the data but if you are in a
hurry, or want to wrap your head around all this to get a better
understanding, here is what you can do :

To be continued… (and where to put some data visualisation + numbers of
line and column for every data.frame to have an idea of what we are
dealing with)
