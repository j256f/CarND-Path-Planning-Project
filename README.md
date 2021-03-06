# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

+On commit 62dc6b I got to code up to what is shown in the Walkthrough video with Aaron Brown and David Silver.
+On commit eadd40 I build on top of what was proposed on the walkthrough to made the first attemp.
+On commit 57a372 I tweak the speed controller algorithm.
+On commit ab4390 I tweaked planner algorithm.

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

#### Navigating Phone Dial


I got this simple idea of dividing the space around the ego vehicle as phone dial sections

![iphone dial](IMG_7625.PNG))


5 being the ego vehicle is at the center and is surrounded by 1 which represents front-left area, 2 representing the area ahead, 3 for the front-right, and so on. 


In this code only 1,2,3,5,7 and 9 areas were used:



#### The Speed controller

The velocity is controlled considering the gap behind a vehicle (Two-sa-gap) that a following vehicle (Five) should keep according to the speed of the vehicle being follwed (Two-sa-v).

So the gap to be minded is the distance that the vehicles Two-sa travels in one second (Two-sa-v = Two-sa-v * 0.4697) 



![speed](HighWay11.jpg)


When ego vehicle (Five) is within a gap distance of the gap, it varies it speed according to the speed difference of the vehicles ((Two-sa-v - Five-v)/Two-sa-v)


#### The Planner, some kind of state machine...

I think half of the job was to generate the arrays with the different vehicles at the different areas around the car. 
First step was to get an array (OneTwoThree) with all vehicles with 300mts ahead of ego vehicle (Five).
Then an array (SevenNine) of all the vehicles within 300 meters behind ego vehicles.
From OneTwoThree three array are generated: One for left lane, Two for lane ahead of Five, and Three, the lane at the right. 
It is important to point out that these lanes are relative to Five.
If Five is at lane#0, Three is goint to be "lane#1 ahead", and Nine be "lane#1 behind".
Two kind of sorted arrays were generated according to linear distance "s" and velocity see lines 801-926 
I am sure there is a better way to generate these arrays, I tried to use Eigen Library but had not enough time to contunue testing.

After the relativistic arrays are generated, I very simple planner can be deduce following three step apprach:
1. What is the best relative lane? add to it the maximum of 1 (only one lane gets maximum, al others 0 at this point) (1109-1011)
2. What the second best? add to it half the maximim (If 1. was not possible, here one lane get the chance to earn some value(113-1120)
3. The third aooportunity is not exclusive, all lanes with vehicles may earn some value.(1127-1129)

After this a pausibility test is taken, if it is not passed, the maximumpenalty is applied (1132-1133)

Then a desambigation test is applied to deduse absolute lane from relative lanes (1163-1182)














