# Synthesis of the Periodic Compensator with Frequency Estimator

<br>

<p align="justify">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In this project, the issue of controlling an industrial robot with two rotary joints in the vertical plane for tracking a periodic reference signal of unknown frequency is elaborated. The kinematic and dynamic model of an industrial robot with two rotary joints in the vertical plane is defined using the Lagrange Euler method. In the case of a constant reference signal, an industrial robot with rotary joints in the vertical plane can be globally asymptotically stabilized using a nonlinear PID controller, but in the case of a periodic reference signal, it is necessary to use a more complex controller. For these needs, a non-linear PID controller with periodic compensator based on the passivity concept is used, which, unlike the periodic compensator based on the internal model concept, has no influence on stability of the closed circuit. The periodic compensator is implemented with a frequency estimator based on the band-stop filter with certain modifications and is analyzed in detail.</p>

<br>


## Project workflow

<br>

<b>Step 1.</b>&nbsp;&nbsp;Designing the model of robot with two rotational joints
<br>
<p align="center"><img src="images%20for%20GitHub/robot%20scheme.jpg" width="320px"></p>
<p align="center"><img src="images%20for%20GitHub/robot%20model.png" width="540px"></p>
<br>

<b>Step 2.</b>&nbsp;&nbsp;Designing the nonlinear PID controller with periodic compensator and frequency estimator
<br>
<p align="center"><img src="images%20for%20GitHub/preiodic%20compensator%20model.png" width="360px"></p>
<br>

<b>Step 3.</b>&nbsp;&nbsp;Implementing the model without controller
<br>
<p align="center"><img src="images%20for%20GitHub/model%20without%20controller.png" width="660px"></p>
<br>

<b>Step 4.</b>&nbsp;&nbsp;Implementing the model with nonlinear PID controller
<br>
<p align="center"><img src="images%20for%20GitHub/model%20with%20controller.png" width="660px"></p>
<br>

<b>Step 5.</b>&nbsp;&nbsp;Implementing the model with nonlinear PID controller and periodic compensator
<br>
<p align="center"><img src="images%20for%20GitHub/model%20with%20controller%20and%20periodic%20compensator.png" width="660px"></p>
<br>

<b>Step 6.</b>&nbsp;&nbsp;Implementing the model with nonlinear PID controller, periodic compensator and frequency estimator
<br>
<p align="center"><img src="images%20for%20GitHub/model%20with%20controller,%20periodic%20compensator%20and%20frequency%20estimator.png" width="660px"></p>
<br>

<b>Step 7.</b>&nbsp;&nbsp;Analyzing the results for different implementations and parameters
<br>
<br>


## Run the project on Windows

<br>

<b>Step 1.</b>&nbsp;&nbsp;Clone the repository:
<pre>
cd %HOMEPATH%

git clone https://github.com/Doc1996/periodic-compensator
</pre>
<br>

<b>Step 2.</b>&nbsp;&nbsp;Run the Matlab script <i>RR_robot_model.m</i> and examine results
<br>
<br>

<b>Step 3.</b>&nbsp;&nbsp;Run the Matlab script <i>RR_robot_RC.m</i> and examine results
<br>
<br>


## Run the project on Linux

<br>

<b>Step 1.</b>&nbsp;&nbsp;Clone the repository:
<pre>
cd $HOME

git clone https://github.com/Doc1996/periodic-compensator
</pre>
<br>

<b>Step 2.</b>&nbsp;&nbsp;Run the Matlab script <i>RR_robot_model.m</i> and examine results
<br>
<br>

<b>Step 3.</b>&nbsp;&nbsp;Run the Matlab script <i>RR_robot_RC.m</i> and examine results
<br>
<br>