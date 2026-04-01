# Method Of Manufactured Solutions
All my work is motivated by this [paper](./resources/mms-paper.pdf)


So we have a basic set up working in the periodic-bc folder. It can be found [here](../MMS/periodic-bc/research-MMS-mbailly.pdf)

I have found this [paper](./resources/sandia-paper.pdf) very interesting. It contains a list of properties the Manufactured Solution should adhere to. It also have some cool examples. Let's try to implement some! 




- Motivation Slide 1
- - Finding Bugs in Numerics can be very challenging.
- - It is helpful to have a true solution to reference.
- Motivation Slide 2 (maybe not needed)
- - I'll insert two graphs here showing how some methods can blow up and it is obvious that they are incorrect while other methods can look correct but still be off.
- Code Verification Slide 1
- - Many different definitions
- - Convincingly demonstrate that the code is correctly solving PDEs and boundary conditions
- Code Verification Slide 1
- - The easiest and most convincing method: know the true solution. Of course, we don't always know the true solution
- - Only tests the correctness of the code, doesn't debug the code for you
- Method of Manufactured Solutions (MMS) Slide 1
- - The user manufactures a solution
- - Solves their desired PDE with an additional forcing term
- Method of Manufactured Solutions (MMS) Slide 2
- - We want to solve the PDE $L(u(t,\mathbf{x})) = 0$
- - Create a manufactured solution $U(t,\mathbf{x})$.
- - Plug the manufactured solution into PDE $L(U(t,\mathbf{x})) = Q(t,\mathbf{x}).$
- - Now $U(t,\mathbf{x})$ is a solution to $L(u(t,\mathbf{x})) = Q(t,x)$ 
- Example
- - put the advection example here.
- Method of Manufactured Solutions (MMS) Slide 3
- - Now knowing the solution $U(t,\mathbf{x})$, verify the correctness of your method on the PDE $L(u(t,\mathbf{x})) = Q(t,x)$
- - Compare numerical result to manufactured solution
- - Verify order of accuracy (will be the same order as the method without the forcing term)
- Guidelines for Creating Manufactured Solution
- - Sufficiently general to exercise all terms in the governing PDE (e.x., if there is time dependence in the PDE, solution also needs time dependence)
- - Meets assumptions for PDE (e.x., if times is expected to be positive, make sure solution is defined for positive)
- - Sufficient number of non-trivial derivatives (e.x., linear solutions won't work for second order PDE)
- - Need not be a physically realistic solution

- Initial and Boundary conditions

- Approximating Accuracy 

- Implementation
- - Add some pseduocode here

- Code and Simple Example

- Future Plans
- - Implement Julia to C code verification
- - Implement different types of boundary conditions
- - Increase complexity of PDEs 
- - - We want to be able to verify Naiver-Stokes!
- - Allow more general solutions
- - - user to define manufactured solution on different domains of 2-D and 3-D space