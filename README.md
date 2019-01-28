## README

all code is in:
```~/analysis/crystal_char```

open up a terminal and `cd` to that directory
```cd ~/analysis/crystal_char```

very soon we will have temperature monitoring.
when it’s ready, you’ll do:
```python auto_process.py —temp```
to enable temp logging while ORCA is taking data for this crystal.

for each crystal,
record serial number (“SNXXXXX”)
take 10 runs  (5 position, 5 voltage)
enter run numbers and serial number into runDB.json:
```atom runDB.json```

then do:
```python auto_process.py -p [SN]```
this creates a folder called [SN] in runDB[“built_path”]
—> make sure you know what I mean by that!

run the Calibration code:
```cd ~/analysis/crystal_char/calibration
./Calibration [path to built directory with SN] [option - pos, volt]
```
For the plots that have to be manually saved,
save them to the built directory that you just created, i.e. runDB[“built_path”] + [SN]

Verify that the files in the built directory are what you expect.
(cd to that directory and look!)

Run the sync to cenpa-rocks and remove both the raw and built files from this machine:
```cd ~/analysis/crystal_char
python auto_process.py -s
```

### for help w/ auto_process, or any other code:
1. actually look at the code!  `atom ~/analysis/crystal_char/auto_process.py`
2. check help: `python auto_process.py -h`
3. google it
4. ask clint

### changing the code

if you make changes to any software: GOOD!  we want your input.
but you have to `git commit` them so that we can track what you’ve done.
this can be done via command line:
```git pull
git add .
git commit -m “very brief message describing my changes”
git push
```
Ideally, every time runDB.json is updated, you should commit.
for now, Clint will probably just do this once a week.*
