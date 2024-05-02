# Empirical bias reducing methods for Structural Equation Models


Last modified: 2024-05-02

## Normal CFA toy example

``` r
# Simulate data
mod <- "fx =~ 1*x1 + 0.8*x2 + 0.6*x3"
n   <- 200
dat <- simulateData(model = mod, sample.nobs = n)
fit <- cfa(model = "fx =~ x1 + x2 + x3", data = dat)
semPaths(fit, "est", rotation = 3)
```

<img src="README_files/figure-commonmark/fig-toyexample-1.png"
id="fig-toyexample" />

Average bias across `B=1000` simulations (sample size n = 200):

    # A tibble: 6 × 5
      par         ml     ebrm      jack     boot
      <chr>    <dbl>    <dbl>     <dbl>    <dbl>
    1 fx=~x2  0.0274  0.00740  0.00664   0.0114 
    2 fx=~x3  0.0133  0.00443  0.00411   0.0168 
    3 x1~~x1 -0.0286  0.00794  0.0144    0.0308 
    4 x2~~x2 -0.0382 -0.0218  -0.0163   -0.0189 
    5 x3~~x3 -0.0134 -0.00554 -0.000257 -0.00865
    6 fx~~fx  0.0181 -0.0185  -0.0150   -0.0255 

## Growth curve models

### Convergence failures and timings

<div>

<div id="gjlbbvafqq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#gjlbbvafqq table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#gjlbbvafqq thead, #gjlbbvafqq tbody, #gjlbbvafqq tfoot, #gjlbbvafqq tr, #gjlbbvafqq td, #gjlbbvafqq th {
  border-style: none;
}
&#10;#gjlbbvafqq p {
  margin: 0;
  padding: 0;
}
&#10;#gjlbbvafqq .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#gjlbbvafqq .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#gjlbbvafqq .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#gjlbbvafqq .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#gjlbbvafqq .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#gjlbbvafqq .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#gjlbbvafqq .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#gjlbbvafqq .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#gjlbbvafqq .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#gjlbbvafqq .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#gjlbbvafqq .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#gjlbbvafqq .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#gjlbbvafqq .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#gjlbbvafqq .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#gjlbbvafqq .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#gjlbbvafqq .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#gjlbbvafqq .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#gjlbbvafqq .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#gjlbbvafqq .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#gjlbbvafqq .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#gjlbbvafqq .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#gjlbbvafqq .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#gjlbbvafqq .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#gjlbbvafqq .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#gjlbbvafqq .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#gjlbbvafqq .gt_left {
  text-align: left;
}
&#10;#gjlbbvafqq .gt_center {
  text-align: center;
}
&#10;#gjlbbvafqq .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#gjlbbvafqq .gt_font_normal {
  font-weight: normal;
}
&#10;#gjlbbvafqq .gt_font_bold {
  font-weight: bold;
}
&#10;#gjlbbvafqq .gt_font_italic {
  font-style: italic;
}
&#10;#gjlbbvafqq .gt_super {
  font-size: 65%;
}
&#10;#gjlbbvafqq .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#gjlbbvafqq .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#gjlbbvafqq .gt_indent_1 {
  text-indent: 5px;
}
&#10;#gjlbbvafqq .gt_indent_2 {
  text-indent: 10px;
}
&#10;#gjlbbvafqq .gt_indent_3 {
  text-indent: 15px;
}
&#10;#gjlbbvafqq .gt_indent_4 {
  text-indent: 20px;
}
&#10;#gjlbbvafqq .gt_indent_5 {
  text-indent: 25px;
}
</style>

| rel | n    | fail_ml | fail_jack | fail_boot | fail_ebrm | time_jack | time_boot | time_ebrm |
|-----|------|---------|-----------|-----------|-----------|-----------|-----------|-----------|
| 0.8 | 15   |         |           | 80.74%    | 0.55%     | 1s        | 25s       | 1s        |
| 0.8 | 20   |         |           | 4.62%     | 0.10%     | 2s        | 52s       | 1s        |
| 0.8 | 50   |         |           |           |           | 5s        | 52s       | 1s        |
| 0.8 | 100  |         |           |           |           | 10s       | 52s       | 1s        |
| 0.8 | 1000 |         |           |           |           | 1m 43s    | 52s       | 1s        |
| 0.5 | 15   |         |           | 80.80%    |           | 1s        | 22s       | 1s        |
| 0.5 | 20   |         |           | 5.47%     |           | 1s        | 43s       | 1s        |
| 0.5 | 50   |         |           |           |           | 4s        | 43s       | 1s        |
| 0.5 | 100  |         |           |           |           | 8s        | 43s       | 1s        |
| 0.5 | 1000 |         |           |           |           | 1m 26s    | 43s       | 1s        |

</div>

</div>

### Relative median bias plot

<img src="README_files/figure-commonmark/fig-bias-growth-1.png"
id="fig-bias-growth" style="width:100.0%" />

## Two factor SEM

TBC
