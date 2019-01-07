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


### Run the program

(Not necessary, data are already retrieved)
```
cd src
python3 alignment.py
python3 space_to_tab.py
```

Boxplots generations:
```
cd src
python3 boxplots.py
```

Statistical tests and heatmaps generation
```
cd src
python3 heatmaps_stats.py
```

## Authors

Master student in bioinformatics at Paris Diderot University.
- [Bonté Pierre-Emmanuel](https://github.com/pebonte/AAindex)
