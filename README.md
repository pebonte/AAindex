# AAindex

## Installation

### Clone the repository
```
git clone https://github.com/pebonte/AAindex.git
```
### Requirements

#### 1. A Linux/MAC OSX distribution.

#### 2. Install the few required Python packages / modules :
```
pip install seaborn
pip install scikit_posthocs
pip install pandas
pip install selenium
pip install matplotlib
pip install numpy
pip install scipy
```
#### 3. A Google chrome webdriver depending on OS.
In the root of the git, a MAC OSX chromedriver is given. You need to download a Linux chromedriver and put it in the root of the git in order to run the program. You can download it [here](https://chromedriver.storage.googleapis.com/index.html?path=2.45/)

#### 4. An Internet connection !!!! (all scripts need internet to run).

### Run the program

#### Convert msf files to plain files
Not necessary, data are already retrieved but you can use it to convert msf files to plain files (easy to parse)
```
cd src
python3 alignment.py
python3 space_to_tab.py
```

#### Boxplots and CSV generations
You can test the script with a few data first by only giving 1 index to work on. An example is given. You just need to manually uncomment this line in the script. 
```
#list_aaindex_ids = ['GRAR740102']
```
To run the script:
```
cd src
python3 boxplots.py
```

#### Statistical tests and heatmaps generation
You can test the script with a few data first by only giving 1 index to work on. An example is given. You just need to manually uncomment this line in the script.
```
#list_aaindex_ids = ['GRAR740102']
```
To run the script:
```
cd src
python3 heatmaps_stats.py
```

## Authors

Master student in bioinformatics at Paris Diderot University.
- [Bonté Pierre-Emmanuel](https://github.com/pebonte/AAindex)
