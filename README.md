# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

On commit XXXX I got to code up to what is shown in the Walkthrough video with Aaron and David Silver.
On commit XXXX I build on top of what was proposed on the walkthrough to made the first attemp
On commit XXXXX I tweak the speed controller algorithm

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

#### Navigating Phone Dial


I got this simple idea of dividing the space around the ego vehicle as phone dial sections

-|-|-
1|2|3
4|5|6
7|8|9

5 being the ego vehicle is at the center and is surrounded by 1 which represents front-left area, 2 representing the area ahead, 3 for the front-right, and so on. 


In this code only the following areas were used:

Here|There|There
-|-|-
1|2|3
|5|
7||9

namely:

Here|There|There
-|-|-
400 mts ahead of left lane|400 mts ahead on same lane as ego vehicle|400 mts ahead of right lane
|ego vehicle|
last 400 mts of left lane||last 400 mts on right line 


#### The Speed controller

The velocity is controlled considering the gap behind a vehicle (Two-sa-gap) that a following vehicle (Five) should keep according to the speed of the vehicle being follwed (Two-sa-v).

So the gap to be minded is the distance that the vehicles Two-sa travels in one second (Two-sa-v = Two-sa-v * 0.4697) 

When ego vehicle (Five) is within a gap distance of the gap, it varies it speed according to the speed difference of the vehicles ((Two-sa-v - Five-v)/Two-sa-v)


#### The Planner










