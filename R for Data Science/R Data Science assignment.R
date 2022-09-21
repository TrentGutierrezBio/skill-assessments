install.packages("tidyverse")
install.packages("plotly")
install.packages("readr")
install.packages("styler")
install.packages(c("colorspace", "corrr",  "cowplot",
                   "ggdark", "ggforce", "ggrepel", "ggridges", "ggsci",
                   "ggtext", "ggthemes", "grid", "gridExtra", "patchwork",
                   "rcartocolor", "scico", "showtext", "shiny",
                   "highcharter", "echarts4r"))
install.packages("kableExtra")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(plotly)
library(styler)
library(kableExtra)


#Problem 1 
data <- (read.csv("gapminder_clean.csv")) %>%
  as.tibble()


data[is.na(data)] <- 0

#Problem 2 


yeardata <- data %>%
  filter(Year == 1962) %>%
  select(Country.Name,CO2.emissions..metric.tons.per.capita., gdpPercap, Year) 


ggplot(yeardata, mapping = aes(x = CO2.emissions..metric.tons.per.capita., y= gdpPercap)) +
  geom_point() 


#Problem 3
yeardata %>% 
  group_by(Year) %>%
  summarize(cor = cor(CO2.emissions..metric.tons.per.capita.,gdpPercap))


p_value <- cor.test(yeardata$CO2.emissions..metric.tons.per.capita., yeardata$gdpPercap)
p_value


#Problem 4
highestdata <- data %>%
  select(Country.Name,Year,gdpPercap,CO2.Emissions = CO2.emissions..metric.tons.per.capita.) %>%
  group_by(Year) %>% 
  summarize(cor = cor(CO2.Emissions, gdpPercap)) %>%
  top_n(1,cor)
highestdata

#Problem 5
scatterplot <- data %>%
  select(Country.Name,Year,CO2.Emission = CO2.emissions..metric.tons.per.capita.,continent,pop,gdpPercap) %>%
  filter(Year == 1962) %>%
  ggplot(mapping = aes(x=CO2.Emission, y=gdpPercap, size = pop, col = continent,)) + 
  geom_point() 
ggplotly(scatterplot)


#########################################################################################################
#Part 2 
#Problem 1 
#What is the relationship between `continent` and `'Energy use (kg of oil equivalent per capita)'`? 

energydata <- data %>%
  select(continent, EnergyUse = Energy.use..kg.of.oil.equivalent.per.capita.) %>%
  group_by(continent) %>%
  summarize(Avg.EnergyUse.in.Kg = mean(EnergyUse))

energydata

workingdata <- c("Africa", "Asia","Americas","Europe","Oceania")

energygraph <- data %>%
  select(Country.Name,Year, continent,EnergyUse = Energy.use..kg.of.oil.equivalent.per.capita.) %>%
  filter(continent %in% workingdata,na.rm = TRUE) %>%
  ggplot(mapping = aes(EnergyUse, col = continent,)) +
  geom_histogram() + 
  scale_x_log10() + 
  facet_wrap(~continent)

ggplotly(energygraph)

#Problem 2
#Is there a significant difference between Europe and Asia with respect to `'Imports of goods and services (% of GDP)'` in the years after 1990? 

FocusContinent <- c("Europe","Asia")
comparedata <- data %>%
  select(continent, Year, ImportEcon = Imports.of.goods.and.services....of.GDP.) %>%
  filter(continent %in% FocusContinent, Year > 1990)

comparedgraph <- ggplot(comparedata, mapping = aes(x=Year, y= ImportEcon), col= continent) +
  geom_smooth() + 
  facet_wrap(~continent)

ggplotly(comparedgraph) 

avgcomparedata <- comparedata %>%
  group_by(continent) %>%
  summarize(average.income = mean(ImportEcon))
avgcomparedata


#3. What is the country (or countries) that has the highest `'Population density (people per sq. km of land area)'` across all years? 
#   (i.e., which country has the highest average ranking in this category across each time point in the dataset?)

#Note: Create a loop to record the country with the highest Population density per year? and then print it out at the end. 
densitydata <- data %>%
  select(Country.Name,Year,Population.density = Population.density..people.per.sq..km.of.land.area.)

year1962data <- densitydata %>%
  filter(Year == 1962) %>% 
  top_n(1 ,Population.density)

year1962data  

X <- 1962 

while(X !=  2012) {
  Neodata <-densitydata %>%
    filter(Year == X) %>%
    top_n(1,Population.density)
    print(paste(Neodata[1,1], "had the record population density of people per sq km of land area of",
                Neodata[1,3], "on the year", Neodata[1,2]))
  X <- X + 5
} 

densitygraph <- ggplot(densitydata, mapping = aes(x=Year,y=Population.density,col= Country.Name)) + 
  geom_point() +
ggplotly(densitygraph)


#4. What country (or countries) has shown the greatest increase in `'Life expectancy at birth, total (years)'` since 1962? 

lifedata <- data %>%
  select(Country.Name,Year, Life.expectancy.yrs = Life.expectancy.at.birth..total..years.) 

toplife <- lifedata %>%
  group_by(Country.Name) %>%
  summarize(Range = max(Life.expectancy.yrs)-min(Life.expectancy.yrs)) %>%
  top_n(5,Range)
toplife

toplifeexp <- c("San Marino","Faroe Islands",
                "Bermuda","Liechtenstein","St. Martin (French part)")

filteredlife <- lifedata %>%
  filter(Country.Name %in% toplifeexp)

lifegraph <- filteredlife %>% 
  ggplot(mapping = aes(x=Year, y= Life.expectancy.yrs, col = Country.Name)) +
  geom_smooth() +
  facet_wrap(~Country.Name)

ggplotly(lifegraph)



################################################################################
gapminder %>%
  kbl() %>%
  kable_styling()













