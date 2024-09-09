# HallThruster.jl version discrepancy observations
Looking for a discrepancy between the most recent version of the code (v0.16.6) and an earlier version of the code (v0.14.3).

## Notes on workflow
To run these tests, I have two Julia environments. Each has Hallthruster.jl installed -- one has the latest version and one has v0.14.3. I boot up julia in this directory and either `]activate newversion` or `]activate oldversion`, depending on which version I'm testing. Then, I `include("test.jl")` to generate the small tables below. I usually run it a second time using the `rl()` function in that file to get a more accurate second timing, neglecting the compilation overhead. This gives me a direct head-to-head way to compare the two versions, with all else held equal.

## Comparison 1 (analytic SPT-100 bfield)

**NOTE:these tests use the analytic spt-100 magnetic field, as I don't have the one in the csv file on this computer**

### v0.16.6
Runtime:           24.35 s
Discharge current: 5.35 A
Ion current:       4.17 A
Thrust:            93.8 mN

### v0.14.3
Runtime:           46.88 s
Discharge current: 5.37 A
Ion current:       3.99 A
Thrust:            93.3 mN

### Observations 
As in my initial tests on this matter, the new version has, if anything, a higher ion current. The discharge current and thrust differ by only O(~1%). The main observable difference is that the new version is nearly twice as fast, which was the reason behind most of the changes made between these two versions.

## Comparison 2 (maagnetic field file)

### v0.16.6

### v0.14.4

### Observations

## Summary
