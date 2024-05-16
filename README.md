# MetaboSpectra

***MetaboSpectra*** is an R package designed for metabolomics spectra match, with some basic functions for mass spectrometry operations built in.

## Input

All ***MetaboSpectra*** operations use standard input in this format:

```R
> data("standardInput")
> standardInput
# A tibble: 2,907 × 6
   feature          precursorMz    rt adduct mz         intensity 
   <chr>                  <dbl> <dbl> <chr>  <list>     <list>    
 1 6.18_170.0924m/z        170.  6.18 [M+H]+ <dbl [18]> <dbl [18]>
 2 6.18_170.0924m/z        170.  6.18 [M+H]+ <dbl [14]> <dbl [14]>
 3 0.47_170.0924m/z        170.  0.47 [M+H]+ <dbl [25]> <dbl [25]>
 4 0.47_170.0924m/z        170.  0.47 [M+H]+ <dbl [20]> <dbl [20]>
 5 1.97_103.0401m/z        103.  1.97 [M-H]- <dbl [5]>  <dbl [5]> 
 6 7.34_301.1798m/z        301.  7.34 [M+H]+ <dbl [86]> <dbl [86]>
 7 2.27_227.0673m/z        227.  2.27 [M-H]- <dbl [15]> <dbl [15]>
 8 2.27_227.0682m/z        227.  2.27 [M-H]- <dbl [55]> <dbl [55]>
 9 2.27_227.0668m/z        227.  2.27 [M-H]- <dbl [12]> <dbl [12]>
10 3.13_228.0979m/z        228.  3.13 [M+H]+ <dbl [30]> <dbl [30]>
# ℹ 2,897 more rows
# ℹ Use `print(n = ...)` to see more rows
```

## Basic operation

The standard input can be converted into a mass spectral matrix (``spMat``)

```R
spMat <- get_spMat(standardRow = standardInput[1, ])
# or
mz <- c(21.3300, 40.1320, 86.3400, 138.3290, 276.5710, 276.5830)
intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
standardRow <- tibble::tibble(mz = list(mz), intensity = list(intensity))
spMat <- get_spMat(standardRow)
```

```R
> spMat
          mz intensity
[1,]  40.132      1300
[2,]  86.340      4030
[3,] 138.329     10000
[4,] 276.571     31600
```

```R
plotSpectra(spMat = spMat)
```

<img src=".\assets\image-20240516091642448.png" alt="image-20240516091642448" style="zoom:67%;" />

Cleaning functions for spectrum, including normalization, merging, and filtering.

```R
spMat <- clean_spMat(spMat)
plotSpectra(spMat = spMat)
```

<img src=".\assets\image-20240516091839588.png" alt="image-20240516091839588" style="zoom:67%;" />

## Spectra match

### Load library

Users can download public libraries and in-house library (TangLab) from ***MetaboLib.ms2*** (Unpublished).

```R
load("hmdbMs2List.RData")
```

### Dot product

```R
searchRes_ndotproduct <- searchLib_ndotproduct(standardInput = standardInput, lib = hmdbMs2List, thread = 8)
```

### Spectral entropy

```R
searchRes_entropy <- searchLib_entropy(standardInput = standardInput, lib = hmdbMs2List, thread = 8)
```

The search library function returns a list whose index corresponds to the input tibble's row index.

```R
> searchRes_entropy[[1]]
# A tibble: 5 × 11
  inchikey     accession adduct collision_energy instrument_type instrument precursorMz predicted mz    intensity score
  <chr>        <chr>     <chr>  <chr>            <chr>           <chr>            <dbl> <lgl>     <lis> <list>    <dbl>
1 JDHILDINMRG… HMDB0000… [M+H]+ 10 V             NA              NA                170. FALSE     <dbl> <dbl [7]> 0.845
2 JDHILDINMRG… HMDB0000… [M+H]+ 20 V             NA              NA                170. FALSE     <dbl> <dbl>     0.841
3 JDHILDINMRG… HMDB0000… [M+H]+ 20 V             NA              NA                170. FALSE     <dbl> <dbl>     0.836
4 JDHILDINMRG… HMDB0000… [M+H]+ 10 V             NA              NA                170. FALSE     <dbl> <dbl [9]> 0.835
5 JDHILDINMRG… HMDB0000… [M+H]+ 35 V             NA              NA                170. FALSE     <dbl> <dbl [7]> 0.833
```

### Results export

```R
idenTibble_entropy <- searchRes2idenTibble(standardInput, searchRes_entropy, top = 5)
```

```R
> idenTibble_entropy
# A tibble: 2,907 × 9
   feature          precursorMz    rt adduct mz         intensity  iden_inchikey                     iden_id iden_score
   <chr>                  <dbl> <dbl> <chr>  <list>     <list>     <chr>                             <chr>   <chr>     
 1 6.18_170.0924m/z        170.  6.18 [M+H]+ <dbl [18]> <dbl [18]> "JDHILDINMRGULE-LURJTMIESA-N"     "HMDB0… "0.844646…
 2 6.18_170.0924m/z        170.  6.18 [M+H]+ <dbl [14]> <dbl [14]> "JDHILDINMRGULE-LURJTMIESA-N"     "HMDB0… "0.862108…
 3 0.47_170.0924m/z        170.  0.47 [M+H]+ <dbl [25]> <dbl [25]> "JDHILDINMRGULE-LURJTMIESA-N"     "HMDB0… "0.850078…
 4 0.47_170.0924m/z        170.  0.47 [M+H]+ <dbl [20]> <dbl [20]> "JDHILDINMRGULE-LURJTMIESA-N"     "HMDB0… "0.838322…
 5 1.97_103.0401m/z        103.  1.97 [M-H]- <dbl [5]>  <dbl [5]>  "SJZRECIVHVDYJC-UHFFFAOYSA-N;AFE… "HMDB0… "0.917271…
 6 7.34_301.1798m/z        301.  7.34 [M+H]+ <dbl [86]> <dbl [86]> ""                                ""      ""        
 7 2.27_227.0673m/z        227.  2.27 [M-H]- <dbl [15]> <dbl [15]> ""                                ""      ""        
 8 2.27_227.0682m/z        227.  2.27 [M-H]- <dbl [55]> <dbl [55]> ""                                ""      ""        
 9 2.27_227.0668m/z        227.  2.27 [M-H]- <dbl [12]> <dbl [12]> "MXHRCPNRJAMMIM-SHYZEUOFSA-N"     "HMDB0… "0.846359…
10 3.13_228.0979m/z        228.  3.13 [M+H]+ <dbl [30]> <dbl [30]> ""                                ""      ""        
# ℹ 2,897 more rows
# ℹ Use `print(n = ...)` to see more rows
```

