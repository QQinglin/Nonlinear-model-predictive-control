# Sensitivity-based distributed optimal control with fixed-point iterations [Project Thesis]

This repository contains implementations of Sensitivity-based distributed optimal control with fixed-point iterations algorithms in Matlab.

## 📌 Project Overview

This project is aimed at understanding and implementing various RL algorithms from scratch, including:

This project implements an algorithm that is a combination of a Fixed-point iteration method
 and a Sensitivity-based algorithm. The algorithm aims to solve the local optimal control
 problem, which is one of the most critically important steps for distributed model predic
tive control. With the result of the algorithm, we can implement the whole distributed
 model predictive control. The combination is numerically efficient in the sense that the
 optimal conditions of the local optimal control problem are solved within a coupled forward
backward integration and the sensitivities can be calculated locally for each neighbor. It
 has one communication step per algorithm iteration. The whole process leads to a fully
 distributed algorithm. Convergence is shown for an upper bound on the horizon time. I
 use the Anderson acceleration to enlarge the original horizon time so that it has better per
formance. And when I implement the distributed model predictive control, the controllers
 have better predictive competence as well. In addition, the algorithm has an important
 trade-off between convergence speed and stability of the sensitivity-based algorithm. And
 except Anderson Acceleration the maximum horizon time length can typically be enlarged
 by damping the iterates as well

## 📁 Folder Structure
  Sensitivity-based distributed optimal control with fixed-point iterations/
├── Bild # results of 5 different models
├── data # results of 5 different models for Latex
├── Model1  # Model with a complex coupling way has 5 agents
├── Model2 # Model with cost function coupling items
├── Model3 # Model with single coupling direction
├── Model4 # Model has nonlinear coupling items in dynamic function
├── References # some referencs in this field
├── Robot # Model of robots
├── README.md






