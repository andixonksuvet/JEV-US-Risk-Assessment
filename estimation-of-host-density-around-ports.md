Density of JEV hosts
================
Dr. Andrea Dixon
2024-07-24

# Background

In order to answer question 50 - What is the probability of infecting a
first local (indigenous) host, given the pathway of entry and the
expected region and time of entry? \[1st transmission step\] - we assume
introduction of the virus will be through an introduced mosquito by air
or inland/coastal port and therefore want to estimate the density of JEV
hosts (feral pigs, domestic pigs, and *ardeidae* birds) in a 7.7 km
radius, the maximum active flight distance for a *Culex* species[^1], of
inland/coastal ports (called seaports in this document) and airports
that had at least one flight/ship from the region of JEV distribution.

# Datasets

``` r
# libraries
library(tidyverse) 
library(sf) 
library(raster) 
library(tigris)
library(knitr)
library(kableExtra)
```

## Location of ports and numbers of ships/flights

See Dixon et al. (doi: TBD) for detailed methods about data sources, but
briefly, flights with an origin in the region of JEV distribution and a
destination in the US were identified from the Bureau of Transportation
(BTS) T-100 Market (all carriers) database for the years 2017-2022[^2].
Airport longitude and latitude was identified as the most extant airport
location from the BTS Master Coordinate file[^3]. Foreign cargo shipping
was identified from US waterway data from the US Army Corp of
Engineers - Waterborne Commerce Statistics Center for the years
2017-2020[^4] (2020 was the most recent year). The number of inbound
ships from a foreign port and containers were collected. Seaport
location was collected from the BTS Principal Port list[^5]. If the port
could not be identified from the list, the location was estimated as the
nearest city.

``` r
airport_mapping <- read_csv("./introduction_map_files/airport_lat_long_merged.csv") %>% 
  mutate(state_abbr = str_split_i(DEST_CITY_NAME, ", ", 2)) # add state abbreviation
seaport_map <- read_csv("./introduction_map_files/shipping_with BTS_2023SEP21_clean_compact.csv")
```

## Feral pigs

Feral pig observations in the US were downloaded from iNaturalist from
Global Biodiversity Information Facility[^6]. Observations were limited
to the time range of 2022 and the U.S and Canada.

``` r
feral_pigs <- read_tsv("./introduction_map_files/0131315-230530130749713.zip") %>% 
  filter(year == 2022) %>% 
  filter(countryCode != "MX")
```

## *Ardeidae* birds

*Ardeidae* bird observations in the US for 2022 were downloaded from
eBird Observation Dataset from Global Biodiversity Information
Facility[^7].

``` r
ardeidae <- read_tsv("./introduction_map_files/0014536-230828120925497.zip") %>% 
  filter(stateProvince != "Midway Islands")
```

# Intersection of ports and hosts

## Ports

``` r
# function to id airports in the different regions at risk
airport.by.region <- function(region_at_risk){
  
  st_as_sf(airport_mapping %>% 
             filter(state_abbr %in% c(region_at_risk)),
           coords = c("LONGITUDE", "LATITUDE"),
           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
}

# function to id seaports to the different regions at risk

seaport.by.region <- function(region_at_risk){
  
  st_as_sf(seaport_map %>% 
             filter(STATE %in% c(region_at_risk)),
           coords = c("Longitude", "Latitude"),
           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
}
```

### Adding index and RAR information

``` r
RAR <-  list(RAR1 = c("ME", "MA", "NH", "NY", "VT"),
             RAR2 = c("CT", "PA", "RI", "IL", "IN", "IA", "KS", 
                      "MI", "MN", "MO", "ND", "NE", "NJ", "OH", 
                      "SD", "WV", "WI"),
             RAR3 = c("AL", "AR", "DE", "FL", "GA", "KY", "LA", 
                      "MD", "MS", "NC", "OK", "SC", "TN", "TX", 
                      "VA", "DC"), # DC added to region 3 (used in feral pig dataset)
             RAR4 = c("AZ", "CO", "MT", "NV", "NM", "UT", "WY"),
             RAR5 = c("CA", "WA", "OR", "ID"),
             RAR6 = c("AK"),
             RAR7 = c("HI")
)

airport_mapping <- airport_mapping %>% 
  mutate(RAR = case_when(
           state_abbr %in% c(RAR$RAR1) ~ "RAR1",
           state_abbr %in% c(RAR$RAR2) ~ "RAR2",
           state_abbr %in% c(RAR$RAR3) ~ "RAR3",
           state_abbr %in% c(RAR$RAR4) ~ "RAR4",
           state_abbr %in% c(RAR$RAR5) ~ "RAR5",
           state_abbr %in% c(RAR$RAR6) ~ "RAR6",
           state_abbr %in% c(RAR$RAR7) ~ "RAR7"),
         
         index = c(1:nrow(airport_mapping))
  )

seaport_map <- seaport_map %>% 
  mutate(RAR = case_when(
           STATE %in% c(RAR$RAR1) ~ "RAR1",
           STATE %in% c(RAR$RAR2) ~ "RAR2",
           STATE %in% c(RAR$RAR3) ~ "RAR3",
           STATE %in% c(RAR$RAR4) ~ "RAR4",
           STATE %in% c(RAR$RAR5) ~ "RAR5",
           STATE %in% c(RAR$RAR6) ~ "RAR6",
           STATE %in% c(RAR$RAR7) ~ "RAR7")
  )
```

## Feral pigs

Observations of feral pigs and port information were converted to `sf`
objects. Each feral pig observation was given a buffer of 600 m[^8],
which corresponds to a core area of 1.09 km<sup>2</sup> for feral
pigs[^9], and determined if it was within 7.7 km of the port location
(see [Background](#background) for rationale).

### Functions and parameters for section

``` r
# function to check distance of ports with 7.7km of the feral pig range (600m)
feralPig.distance.check <- function(port_list, port_type){
  
    switch(port_type,
           
           airport = st_is_within_distance(port_list, 
                        st_buffer(feral_sf, dist = 600),
                        dist = 7700) %>%
             
             set_names(nm = as.character(port_list$DEST_AIRPORT_ID)) %>% 
             
             compact(),
           
           seaport = st_is_within_distance(port_list, 
                        st_buffer(feral_sf, dist = 600),
                        dist = 7700) %>%
             
             set_names(nm = as.character(port_list$index)) %>% 
             
             compact())
                       
}
```

### Summarize

``` r
# putting feral_pigs as a sf

feral_sf <- st_as_sf(feral_pigs, 
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
```

#### Airports

``` r
# create list of airports at sf
airport_lists_RAR <- purrr::map(RAR, ~airport.by.region(region_at_risk = .x))


# airports within 7.7km of the feral pig range (1.3km)?

airports_near_feralPigs <- purrr::map(airport_lists_RAR, ~feralPig.distance.check(port_list = .x, port_type = "airport"))
    


# count of pigs by airport

feralPig_counts <- purrr::map_depth(airports_near_feralPigs, 2,  length) %>% 
  map_df(~list_flatten(.x)) %>% 
  pivot_longer(cols = everything(),
               names_to = "DEST_AIRPORT_ID",
               values_to = "feralPigCount",
               values_drop_na = TRUE) %>% 
  right_join(airport_mapping %>% 
              mutate(DEST_AIRPORT_ID = as.character(DEST_AIRPORT_ID))) %>% 
  mutate(feralPigCount = replace_na(feralPigCount, 0))
```

#### Seaports

``` r
# create list of seaports at sf
seaport_lists_RAR <- purrr::map(RAR, ~seaport.by.region(region_at_risk = .x)) 
```

    ## Warning in min(cc[[1]], na.rm = TRUE): no non-missing arguments to min;
    ## returning Inf

    ## Warning in min(cc[[2]], na.rm = TRUE): no non-missing arguments to min;
    ## returning Inf

    ## Warning in max(cc[[1]], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in max(cc[[2]], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

``` r
# seaports within 7.7km of the feral pig range (1.3km)?

seaports_near_feralPigs <- purrr::map(seaport_lists_RAR, ~feralPig.distance.check(port_list = .x, port_type = "seaport"))

# count of pigs by seaport

feralPig_counts_sea <- purrr::map_depth(seaports_near_feralPigs, 2,  length) %>% 
  map_df(~list_flatten(.x)) %>% 
  pivot_longer(cols = everything(),
               names_to = "index",
               values_to = "feralPigCount",
               values_drop_na = TRUE) %>% 
  right_join(seaport_map %>% 
               mutate(index = as.character(index))) %>% 
  mutate(feralPigCount = replace_na(feralPigCount, 0))
```

## Domestic pigs

As we do not have exact location of pig farms, the density of pig
livestock in the county of the port of interest will be used as a proxy
for domestic pig hosts in proximity to the site of introduction. Density
was calculated from hog production from the NASS 2017 Ag census[^10] and
the square kilometer of each county based on the tigris package[^11].

### Summarize

``` r
# read in file from "plotting domestic pig density" script
domestic_hog <- readRDS("./introduction_map_files/domestic_hog_density_county.rds") %>% 
  st_transform(crs = "WGS84") #NADS83 to WGS84 to match airport mapping
```

#### Airports

``` r
# retrieving the county where each airport is in
airport_domestic <- st_contains(domestic_hog, 
                               st_as_sf(airport_mapping,
                                        coords = c("LONGITUDE", "LATITUDE"),
                                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
                               ) %>% 
  
  set_names(nm = as.character(domestic_hog$COUNTYNS)) %>% 
  
  compact() 


# creating dataframes and binding into a df
airport_domestic_df <- map(airport_domestic, 
                           ~tibble(COUNTYNS = names(.x),
                                   airport_index = pluck(.x))
                           ) %>% 
  
  bind_rows(.id = "COUNTYNS") %>%
  
  left_join(domestic_hog %>% 
              dplyr::select(c(COUNTYNS,
                              STATEFP,
                              No.Hog,
                              hog_production,
                              hog_density,
                              hog_binned)
                            ),
            by = "COUNTYNS") %>%
  
  left_join(airport_mapping, by = c("airport_index" = "index")
            ) %>%
  
  mutate(hog_binned = fct_na_value_to_level(hog_binned, level = "None"))
```

#### Seaports

``` r
# retrieving the county where each seaport is in
seaport_domestic <- st_contains(domestic_hog, 
                               st_as_sf(seaport_map,
                                        coords = c("Longitude", "Latitude"),
                                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
                               ) %>% 
  
  set_names(nm = as.character(domestic_hog$COUNTYNS)) %>% 
  
  compact() 


# creating dataframes and binding into a df
seaport_domestic_df <- map(seaport_domestic, 
                           ~tibble(COUNTYNS = names(.x),
                                   seaport_row_number = pluck(.x))
                           ) %>% 
  
  bind_rows(.id = "COUNTYNS") %>%
  
  left_join(domestic_hog %>% 
              dplyr::select(c(COUNTYNS,
                              STATEFP,
                              No.Hog,
                              hog_production,
                              hog_density,
                              hog_binned)
                            ),
            by = "COUNTYNS") %>%
  
  left_join(seaport_map %>% 
              mutate(row_number = c(1:nrow(seaport_map))), 
            by = c("seaport_row_number" = "row_number")
            ) %>%
  
  mutate(hog_binned = fct_na_value_to_level(hog_binned, level = "None"))
```

## *Ardeidae* birds

### Functions and parameters for section

``` r
# RAR list with surrounding states as necessary from visual assessment
bird_RAR <-  list(RAR1_B = c("Maine", "Massachusetts", "New Hampshire", "New York", "Vermont",
                           "Connecticut", "Rhode Island", "New Jersey", "Pennsylvania"),
                  RAR2_B = c("Massachusetts", "New Hampshire", "New York", "Vermont",
                           "Connecticut", "Rhode Island", "New Jersey", "Pennsylvania",
                           "Illinois", "Indiana", "Iowa", "Kansas", "Michigan", 
                           "Minnesota", "Missouri", "North Dakota", "Nebraska", "New Jersey", 
                           "Ohio", "South Dakota", "West Virginia", "Wisconsin", "Kentucky", 
                           "Oklahoma"),
                  RAR3_B = c("Alabama", "Arkansas", "Delaware", "Florida", "Georgia", "Kentucky",
                           "Louisiana", "Maryland", "Mississippi", "North Carolina", "Oklahoma", 
                           "South Carolina", "Tennessee", "Texas", "Virginia", 
                           "District of Columbia", "Indiana", "New Mexico"), 
                  RAR4_B = c("Arizona", "Colorado", "Montana", "Nevada", "New Mexico", "Utah", "Wyoming", "California"),
                  RAR5_B = c("California", "Washington", "Oregon", "Idaho"),
                  RAR6_B = c("Alaska"),
                  RAR7_B = c("Hawaii")
)

# Months conducive/optimal for mosquito based on Ana's notes
RAR_vector_season <- list(RAR1_v = c(5,6,7,8),
                          RAR2_v = c(4,5,6,7,8,9),
                          RAR3_v = c(2,3,4,5,6,7,8,9,10),
                          RAR4_v = c(3,4,5,6,7,8,9,10),
                          RAR5_v = c(3,4,5,6,7,8,9,10),
                          RAR6_v = c(6,7,8),
                          RAR7_v = c(1,2,3,4,5,6,7,8,9,10,11,12)
                          )

# function to check distance of ports with 7.7km of the bird observations
ardeidae.distance.check <- function(port_list, bird_states, port_type){
  
  switch(port_type,
         
         airport = st_is_within_distance(port_list, bird_states, dist = 7700) %>% 
           
           set_names(nm = as.character(port_list$DEST_AIRPORT_ID)) %>% 
           
           compact(),
         
         seaport = st_is_within_distance(port_list, bird_states, dist = 7700) %>% 
           
           set_names(nm = as.character(port_list$index)) %>% 
           
           compact()
         
  )
  
}


# function to id bird in the different regions at risk 
bird.by.region <- function(region_at_risk){
  
  st_as_sf(ardeidae %>% 
             filter(stateProvince %in% c(region_at_risk)),
           coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
}
```

### Summarize

``` r
# bird obs broken up into regions for speed
bird_region <- purrr::map(bird_RAR, ~ bird.by.region(region_at_risk = .x))
```

#### Airports

``` r
# airports within 7.7km of an ardeidae observation?

airports_near_ardeidae <- purrr::map2(airport_lists_RAR, bird_region, 
                                      ~ ardeidae.distance.check(port_list = .x, 
                                                                bird_states = .y,
                                                                port_type = "airport")
                                      ) 



# count of birds by airport

ardeidae_counts <- purrr::map_depth(airports_near_ardeidae, 2,  length) %>% 
  map_df(~list_flatten(.x)) %>% 
  pivot_longer(cols = everything(),
               names_to = "DEST_AIRPORT_ID",
               values_to = "ardeidaeCount",
               values_drop_na = TRUE) %>% 
  right_join(airport_mapping %>% 
              mutate(DEST_AIRPORT_ID = as.character(DEST_AIRPORT_ID))
             ) %>% 
  mutate(ardeidaeCount = replace_na(ardeidaeCount, 0))
```

#### Seaports

``` r
# seaports within 7.7km of an ardeidae observation?

seaports_near_ardeidae <- purrr::map2(seaport_lists_RAR, bird_region, 
                                      ~ ardeidae.distance.check(port_list = .x, 
                                                                bird_states = .y,
                                                                port_type = "seaport")
                                      ) 



# count of birds by seaport

ardeidae_counts_sea <- purrr::map_depth(seaports_near_ardeidae, 2,  length) %>% 
  map_df(~list_flatten(.x)) %>% 
  pivot_longer(cols = everything(),
               names_to = "index",
               values_to = "ardeidaeCount",
               values_drop_na = TRUE) %>% 
  right_join(seaport_map %>% 
              mutate(index = as.character(index))
             ) %>% 
  mutate(ardeidaeCount = replace_na(ardeidaeCount, 0))
```

## Summary

### Airports

Information on airports in “Non-disclosed” counties:

- Washington (53): Moses Lake - 20,000 animals in county (based on
  <https://www.nass.usda.gov/Publications/AgCensus/2017/Online_Resources/Ag_Atlas_Maps/17-M211g.php>)

- Ohio (39): Columbus - at edge of heavy production

- Kansas (20): Wichita - amongst counties with moderate production

``` r
airport_domestic_df %>% 
  group_by(RAR) %>% 
  count(hog_binned) %>% 
  pivot_wider(names_from = "hog_binned",
              values_from = "n",
              values_fill = 0) %>% 
  dplyr::select(c(RAR, None, `< 1 - 10`, `11 - 50`, `Not disclosed`)) %>% 
  kable(caption = "The number of airports by hog county density and RAR") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
The number of airports by hog county density and RAR
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
None
</th>
<th style="text-align:right;">
\< 1 - 10
</th>
<th style="text-align:right;">
11 - 50
</th>
<th style="text-align:right;">
Not disclosed
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

``` r
feralPig_counts %>% 
  group_by(RAR) %>% 
  summarise("Number of airports" = n_distinct(DEST_AIRPORT_ID),
            "Mean obs Feral pig" = mean(feralPigCount),
            "SD" = sd(feralPigCount)/sqrt(n()),
            "Median obs feral pig" = median(feralPigCount),
            "range obs feral pig" = paste0(min(feralPigCount), " - ", max(feralPigCount))
            ) %>% 
kable(digits = 2,
      caption = "Summary table of the number of feral pig observations within 7.7 km of an international airport with at least one plane coming from the JEV region of distribution") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
Summary table of the number of feral pig observations within 7.7 km of
an international airport with at least one plane coming from the JEV
region of distribution
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
Number of airports
</th>
<th style="text-align:right;">
Mean obs Feral pig
</th>
<th style="text-align:right;">
SD
</th>
<th style="text-align:right;">
Median obs feral pig
</th>
<th style="text-align:left;">
range obs feral pig
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 3
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR4
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
</tbody>
</table>

``` r
ardeidae_counts %>% 
  group_by(RAR) %>% 
  summarise("Number of airports" = n_distinct(DEST_AIRPORT_ID),
            "Mean obs ardeidae" = mean(ardeidaeCount),
            "SEM" = sd(ardeidaeCount)/sqrt(n()),
            "Median obs ardeidae" = median(ardeidaeCount),
            "range obs ardeidae" = paste0(min(ardeidaeCount), " - ", max(ardeidaeCount))
            ) %>% 
kable(digits = 2,
      caption = "Summary table of the number of Ardeidae observations within 7.7 km of an international airport with at least one plane coming from the JEV region of distribution") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
Summary table of the number of Ardeidae observations within 7.7 km of an
international airport with at least one plane coming from the JEV region
of distribution
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
Number of airports
</th>
<th style="text-align:right;">
Mean obs ardeidae
</th>
<th style="text-align:right;">
SEM
</th>
<th style="text-align:right;">
Median obs ardeidae
</th>
<th style="text-align:left;">
range obs ardeidae
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1072.00
</td>
<td style="text-align:right;">
467.64
</td>
<td style="text-align:right;">
319.0
</td>
<td style="text-align:left;">
102 - 4645
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
390.24
</td>
<td style="text-align:right;">
174.19
</td>
<td style="text-align:right;">
80.0
</td>
<td style="text-align:left;">
1 - 4377
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
442.00
</td>
<td style="text-align:right;">
130.97
</td>
<td style="text-align:right;">
292.5
</td>
<td style="text-align:left;">
5 - 2863
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR4
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
286.56
</td>
<td style="text-align:right;">
127.16
</td>
<td style="text-align:right;">
126.0
</td>
<td style="text-align:left;">
2 - 1186
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
1963.81
</td>
<td style="text-align:right;">
624.32
</td>
<td style="text-align:right;">
568.0
</td>
<td style="text-align:left;">
1 - 9363
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:left;">
0 - 2
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
217.00
</td>
<td style="text-align:right;">
86.23
</td>
<td style="text-align:right;">
221.0
</td>
<td style="text-align:left;">
4 - 422
</td>
</tr>
</tbody>
</table>

### Seaports

- Florida - in a county w. 20K hogs according to APHIS map
- Pennsylvania - border of high density areas
- Virginia - border of higher density areas

``` r
seaport_domestic_df %>% 
  group_by(RAR) %>% 
  count(hog_binned) %>% 
  pivot_wider(names_from = "hog_binned",
              values_from = "n",
              values_fill = 0) %>% 
  dplyr::select(c(RAR, None, `< 1 - 10`, `11 - 50`, `Not disclosed`)) %>% 
  kable(caption = "The number of seaports by hog county density and RAR") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
The number of seaports by hog county density and RAR
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
None
</th>
<th style="text-align:right;">
\< 1 - 10
</th>
<th style="text-align:right;">
11 - 50
</th>
<th style="text-align:right;">
Not disclosed
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
</tr>
</tbody>
</table>

``` r
feralPig_counts_sea %>% 
  group_by(RAR) %>% 
  summarise("Number of airports" = n_distinct(index),
            "Mean obs Feral pig" = mean(feralPigCount),
            "SEM" = sd(feralPigCount)/sqrt(n()),
            "Median obs feral pig" = median(feralPigCount),
            "range obs feral pig" = paste0(min(feralPigCount), " - ", max(feralPigCount))
            ) %>% 
kable(digits = 2,
      caption = "Summary table of the number of feral pig observations within 7.7 km of a seaport with at least one cargo ship with an origin port in the JEV region of distribution") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
Summary table of the number of feral pig observations within 7.7 km of a
seaport with at least one cargo ship with an origin port in the JEV
region of distribution
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
Number of airports
</th>
<th style="text-align:right;">
Mean obs Feral pig
</th>
<th style="text-align:right;">
SEM
</th>
<th style="text-align:right;">
Median obs feral pig
</th>
<th style="text-align:left;">
range obs feral pig
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 1
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 0
</td>
</tr>
</tbody>
</table>

``` r
ardeidae_counts_sea %>% 
  group_by(RAR) %>% 
  summarise("Number of airports" = n_distinct(index),
            "Mean obs ardeidae" = mean(ardeidaeCount),
            "SEM" = sd(ardeidaeCount)/sqrt(n()),
            "Median obs ardeidae" = median(ardeidaeCount),
            "range obs ardeidae" = paste0(min(ardeidaeCount), " - ", max(ardeidaeCount))
            ) %>% 
kable(digits = 2,
      caption = "Summary table of the number of Ardeidae observations within 7.7 km of a seaport with at least one cargo ship with an origin port in the JEV region of distribution") %>% 
  kable_paper()
```

<table class=" lightable-paper" style="color: black; font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;">
<caption>
Summary table of the number of Ardeidae observations within 7.7 km of a
seaport with at least one cargo ship with an origin port in the JEV
region of distribution
</caption>
<thead>
<tr>
<th style="text-align:left;">
RAR
</th>
<th style="text-align:right;">
Number of airports
</th>
<th style="text-align:right;">
Mean obs ardeidae
</th>
<th style="text-align:right;">
SEM
</th>
<th style="text-align:right;">
Median obs ardeidae
</th>
<th style="text-align:left;">
range obs ardeidae
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RAR1
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
1224.91
</td>
<td style="text-align:right;">
475.80
</td>
<td style="text-align:right;">
903
</td>
<td style="text-align:left;">
5 - 5145
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR2
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
1173.88
</td>
<td style="text-align:right;">
310.96
</td>
<td style="text-align:right;">
570
</td>
<td style="text-align:left;">
24 - 5145
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR3
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
707.26
</td>
<td style="text-align:right;">
109.81
</td>
<td style="text-align:right;">
477
</td>
<td style="text-align:left;">
0 - 2982
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR5
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
1252.90
</td>
<td style="text-align:right;">
222.89
</td>
<td style="text-align:right;">
527
</td>
<td style="text-align:left;">
86 - 3893
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR6
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
3.29
</td>
<td style="text-align:right;">
1.85
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
0 - 28
</td>
</tr>
<tr>
<td style="text-align:left;">
RAR7
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
303.40
</td>
<td style="text-align:right;">
89.93
</td>
<td style="text-align:right;">
273
</td>
<td style="text-align:left;">
85 - 630
</td>
</tr>
</tbody>
</table>

``` r
sessionInfo()
```

    ## R version 4.3.3 (2024-02-29 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 11 x64 (build 22621)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] kableExtra_1.4.0 knitr_1.46       tigris_2.1       raster_3.6-26   
    ##  [5] sp_2.1-3         sf_1.0-16        lubridate_1.9.3  forcats_1.0.0   
    ##  [9] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.2      readr_2.1.5     
    ## [13] tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.0    tidyverse_2.0.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.4       xfun_0.43          lattice_0.22-5     tzdb_0.4.0        
    ##  [5] vctrs_0.6.5        tools_4.3.3        generics_0.1.3     parallel_4.3.3    
    ##  [9] proxy_0.4-27       fansi_1.0.6        highr_0.10         pkgconfig_2.0.3   
    ## [13] KernSmooth_2.23-22 uuid_1.2-0         lifecycle_1.0.4    compiler_4.3.3    
    ## [17] munsell_0.5.1      terra_1.7-71       codetools_0.2-19   htmltools_0.5.8.1 
    ## [21] class_7.3-22       yaml_2.3.8         crayon_1.5.2       pillar_1.9.0      
    ## [25] classInt_0.4-10    wk_0.9.1           tidyselect_1.2.1   digest_0.6.35     
    ## [29] stringi_1.8.3      fastmap_1.1.1      grid_4.3.3         colorspace_2.1-0  
    ## [33] cli_3.6.2          magrittr_2.0.3     utf8_1.2.4         e1071_1.7-14      
    ## [37] withr_3.0.0        scales_1.3.0       rappdirs_0.3.3     bit64_4.0.5       
    ## [41] timechange_0.3.0   rmarkdown_2.26     httr_1.4.7         bit_4.0.5         
    ## [45] hms_1.1.3          evaluate_0.23      viridisLite_0.4.2  s2_1.1.6          
    ## [49] rlang_1.1.3        Rcpp_1.0.12        glue_1.7.0         DBI_1.2.2         
    ## [53] xml2_1.3.6         svglite_2.1.3      rstudioapi_0.16.0  vroom_1.6.5       
    ## [57] R6_2.5.1           systemfonts_1.0.6  units_0.8-5

[^1]: Verdonschot, PFM and Besse-Lototskaya, AA (2014) Flight distance
    of mosquitoes (Culicidae): A metadata analysis to support the
    management of barrier zones around rewetted and newly constructed
    wetlands. Limnologica, 45: 69-79. doi:
    <https://doi.org/10.1016/j.limno.2013.11.002>

[^2]: <https://transtats.bts.gov/TableInfo.asp?gnoyr_VQ=FMF&QO_fu146_anzr=Nv4%20Pn44vr45&V0s1_b0yB=D>;
    accessed Aug 22, 2022

[^3]: <a
    href="https://www.transtathttps://www.transtats.bts.gov/Fields.asp?gnoyr_VQ=FLLs.bts.gov/tables.asp?gnoyr_VQ=FLL&amp;flf_gnoyr_anzr=g_ZNfgRe_PbeQ"
    class="uri">https://www.transtats.bts.gov/Fields.asp?gnoyr_VQ=FLL</a>;

[^4]: <https://publibrary.planusace.us/#/series/Waterborne%20Foreign%20Cargo>;
    accessed Aug 25, 2022

[^5]: <https://geodata.bts.gov/datasets/usdot>[::principal-ports/about](https://geodata.bts.gov/datasets/usdot::principal-ports/about)
    ; accessed Sep 19, 2023

[^6]: GBIF.org (2 August 2023) GBIF Occurrence Download
    <https://doi.org/10.15468/dl.aydh84>

[^7]: GBIF.org (12 September 2023) GBIF Occurrence Download
    <https://doi.org/10.15468/dl.g8dmec>

[^8]: This buffer was included as pigs are not stationary. However, the
    removal of this buffer does not change the results of this analysis
    (results not shown). As the buffer size increases, theoretically,
    more pigs would be within the 7.7km range of a port.

[^9]: Gabor, TM, Hellgren, EC, Bussche, RA, Silvy, NJ. (1999)
    Demography, sociospatial behaviour and genetics of feral pigs (Sus
    scrofa) ina semi-arid environment. Journal of Zoology 247(3):
    311-322.

    Franckowiak and Poché (2018;
    (<https://doi.org/10.1674/0003-0031-179.1.28>) also report kernal
    density estimates of feral hog core area sizes of 1.16 (+/- 0.95)
    and 1.51 (+/- 0.96) km<sup>2</sup> for two field sites in Texas.

[^10]: <https://www.nass.usda.gov/Publications/AgCensus/2017/index.php>

[^11]: Walker K (2023). *tigris: Load Census TIGER/Line Shapefiles*. R
    package version 2.0.1, <https://CRAN.R-project.org/package=tigris>.
