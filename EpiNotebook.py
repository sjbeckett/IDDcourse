import marimo

__generated_with = "0.10.8-dev3"
app = marimo.App()


@app.cell(hide_code=True)
def _():
    import marimo as mo #load the marimo package

    #generate some text output in markdown
    mo.md(
        f"""
        # Epidemiological Modeling Exercises
        #### Dr. Stephen Beckett, University of Maryland (beckett@umd.edu)
        # 
        **Welcome** to the simulation exercise section for Day 2 of BSCI439C /    BIOL708F: "Infectious disease dynamics: a systems approach". Well done on opening the simulation notebook! 🙌 
        In this notebook are a series of interactive exercises and code to explore some  epidemiological models. We are using a notebook approach as we want to show you the codes that such models are built upon -- but in a more user friendly environment, with legible text and graphics, and interactive elements, which we hope will be fun! The key point here is not in learning how to code, which is challenging and goes beyond what we offer in this course, but to see how concepts we have discussed can be translated into models and simulations. Here's a test block -- click play (▶) near the bottom right of the cell below to calculate 5*7. 
        """
        )
    return (mo,)


@app.cell
def _():
    5*7
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        We will start with the simple SIR model that we talked about in class, and work through some exercises that build upon this framework to add some additional concepts. You should gain an appreciation for how computational simulations of epidemics are built, what assumptions are built into them, and how changes to the underlying model parameters and model structures can shape epidemic trajectories. As you work through - we ask you to reflect. Consider asking questions such as: what are the assumptions that each model is making, and how realistic might they be? how sensitive are model outcomes to the model parameters?

        Do not hesitate to ask Stephen and Gabi questions, we do not want you to get stuck. 🧗🚧🙋

        There are 5 different exercises to work through in this interactive notebook. Have fun!

        ## About this notebook
        This notebook was designed and written using a marimo.io notebook which is coded in python - a general computing programming language. The bonus of marimo is that it can be used as an interactive notebook environment with reactivity! This means that we can use objects such as sliders to change input parameters, that will automatically update outputs - such as figures showing simulation data. This makes the notebook easy to use -- even without needing to learn computer programming!

        [ A sidenote: The original exercises were coded in Julia, using the Pluto.jl library in the Julia computer programming language before being ported to python (marimo was inspired by Pluto.jl amongst other tools). Julia is nearly 13 years old (much younger than other languages you may have heard of such as python, or R!) and was designed with scientific computing in mind. Stephen uses Julia in his research. You can read more about it here: [https://julialang.org/](https://julialang.org/) ]

        ## Refreshers:
        **State variables** - represent the quantities that we want to track changes of through time e.g. S, I and R in the SIR model.

        **Parameters** - represent quantities that we assume are constant and act to modify state variables e.g., transmission rate β in the SIR model.

        **Differential Equation** - A set of equations that describe how quickly state variables change with respect to another quantity (here, we will always consider changes with respect to time).

        # 1) The SIR model

        The **S**usceptible-**I**nfected-**R**emoved (SIR) model, is perhaps the most  well known epidemiological model. The model considers that individuals in a population can be categorized as belonging to one of three states: those who are susceptible to being infected (S), those who are infectious (I), and those who are removed from the infection chain (R). The removed category includes those who have recovered from infection and have gained immunity against infection (how it is classically considered - and how we will generally consider it here), but could also include those who are resistant to infecton, or those whose infections led to mortality. The model supposes that those who are infectious can spread the disease to those who are susceptible with a transmission rate given by β; and that those who are infectious recover at a rate given by γ. This model can be written as a set of coupled differential equations to track the changes in the population fractions S, I, and R through time as:

        $$\frac{dS}{dt} = \overset{transmission}{\overbrace{-βSI}}$$

        $$\frac{dI}{dt} = \overset{transmission}{\overbrace{βSI}} - \overset{recovery}{\overbrace{γI}}$$

        $$\frac{dR}{dt} = \overset{recovery}{\overbrace{γI}}.$$

        **Note 1:** this model shows a flow from S through I through R.
        **Note 2:** As we are tracking population fractions, note that, for all times t we expect that $S(t) + I(t) + R(t) =1$. This condition can be a nice diagnostic to check this model works as intended.

        Now let's look at how we might turn this mathematical model into computer code and simulation. First we need to load in some packages and libraries that contain codes that will help us run the simulations, plot the outputs, and build the interactivity of this notebook:
        """
    )
    return


@app.cell
def _():
    #import some packages with functions to help run the notebook
    import matplotlib.pyplot as plt   #helps with plotting
    import numpy as np    #helps with arrays
    from scipy.integrate import odeint   #helps with numerical integration of dynamical systems
    return np, odeint, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        First, we will define the SIR model as a function which is designed to be compatible with the odeint function in the scipy package - which is used to simulate differential equations. In the below code block:

        * the rate of change of the state variables is defined by du
        * the state variables are represented by u
        * the parameters by p, and
        * time (which we do use here) is given by t.
        """
    )
    return


@app.cell
def _():
    # The SIR model differential equations.
    def SIRmodel(u, t, beta, gamma):
        S, I, R = u
        dSdt = -beta * S * I
        dIdt = beta * S * I - gamma * I
        dRdt = gamma * I
        return dSdt, dIdt, dRdt
    return (SIRmodel,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""Ok! Now we have our model — we now need to consider a particular scenario to numerically simulate and write out the corresponding parameters. In the below we specify a simulation over 100 days, where 1 in 10,000 people in the population are initially infected by a disease with transmission rate of β = 0.4/day and remains infectious for an average of 10 days i.e., the recovery rate γ is 1/10 = 0.1/day.""")
    return


@app.cell
def _(SIRmodel, beta, gamma, np, odeint):
    # Total population, N.
    N = 10000
    # Initially infected and recovered individuals at time 0, I0 and R0.
    I0, R0 = 1/N, 0
    # Initial susceptible population at time 0.
    S0 = 1 - I0 - R0

    # Initial conditions vector
    u0 = S0, I0, R0

    #simulation timing vector
    t = np.linspace(0, 200, 200) # creates a vector from 0 to 200, with 200 elements

    # Integrate the SIR equations over the time grid, t.
    ret = odeint(SIRmodel, u0, t, args=(beta.value, gamma.value))
    S, I, R = ret.T
    return I, I0, N, R, R0, S, S0, ret, t, u0


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        #Almost there!
        Now we need to consider the properties of the rates of transmission (β) and recovery (γ), which will differ depending on the disease. The values of these important parameters can really change the expected disease dynamics. Move the sliders below, and see how it changes the disease dynamics in the figure below the sliders!

    	What happens when transmission rates are fast vs. slow?              
    	What happens when recovery rates are fast vs. slow?
        """
    ).callout("info")
    return


@app.cell(hide_code=True)
def _(mo):
    beta = mo.ui.slider(0.01, 2, value=0.4, step=0.01, label='β')
    gamma = mo.ui.slider(0.01, 2, value=0.1, step=0.01, label='γ')
    mo.md(
        f"""
        **Epidemiological parameters.**

        {beta}
        {gamma}
        """
    )
    return beta, gamma


@app.cell
def _(I, R, S, plt, t):
    plt.figure(figsize=(8,4))
    plt.plot(t,S,label='S')
    plt.plot(t,I,label='I')
    plt.plot(t,R,label='R')
    plt.xlabel('Time (days)')
    plt.ylabel('Population fraction')
    plt.legend()
    plt.title("Simulation of epidemiological population dynamics with the SIR model")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ### Question:
        How many different categories of dynamic outcomes are there in the SIR model?
        """
    ).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        With some of the above parameter combinations the disease can take off and many (sometimes nearly all) people can become infected. But, this is not always the case. For some combinations of parameters the disease fails to spread much further than the intial infected population fraction. The SIR model can represent an epidemic, or a failed epidemic. Something it cannot represent is a scenario of a disease that becomes sustained in a population i.e., an endemic disease. We will return to this later.

        In fact, there is a metric for understanding whether an initial outbreak is likely to spread in a population that we can calculate. This is known as the **basic reproduction number, $\cal{R}_0$,** which represents the average number of new infections likely to be generated by one infectious individual in an otherwise susceptible population. If the basic reproduction number is less than one it means that, on average, each new infection will lead to fewer infections in the future - the disease is unlikely to lead to an outbreak. On the other hand, if the basic reproduction number is greater than one, each new infection will generate even more infections and an outbreak is likely.

        The basic reproduction number is calculated by multiplying the rate of transmission by infectious individuals by the average duration that they are infectious (remember the inverse of a rate is a duration). In the SIR model that is:

        $$\mathcal{R}_0 = \frac{\beta}{\gamma}.$$

        Go back to the simulation above -- how does the simulation behave as you change the parameters relative to $\mathcal{R}_0 = 1$ ?

        A related quantity of interest is the **effective reproduction number, $\mathcal{R}_e$**, that accounts for the changes in underlying susceptibility of the population. As the number of susceptibles decrease there are fewer opportunities for an infection to spread, as infected individuals will be more likely to interact with individuals who are also infected, or recovered. This explains why it is *only* when β is much greater than γ that we observe populations in which most people have been infected. This effect can be accounted for if we track the fraction of remaining susceptibles (who in the future could become infectious) through time as:

        $$\mathcal{R}_e(t) = \frac{\beta}{\gamma}.S(t).$$

        As the susceptible population depletes, the effective reproduction number will decrease to less than one, and the number of infected individuals will decline. This represents a concept termed **herd immunity**, which is the level of immunity a population needs to acquire (via infection, or vaccination - though this is not included in the SIR model) in order to stop an outbreak from spreading.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        ### Question:
        If $S+I+R=1$, we know $R$ are immune, and we assume that $I$ is approximately $0$, then what fraction of the population needs to be in $R$ such that the herd immunity threshold is reached?"""
    ).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # 2) Plot interactivity
        Did you notice the interactive elements of the outbreak simulation above? If you go back to the plot above and click the three ... button near the play button, you will open a menu on actions - including one that allows you to download a figure ("Export output as PNG")!

        Go back and try this functionality!
        Note - this may not work in 'app mode'.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # 3) Endemic outbreaks and waning immunity
        As mentioned previously, the SIR model can not be used to represent diseases that are persistent (i.e., those in which the number of infectious individuals stays greater than 0). This is due to the depletion of susceptible individuals in the SIR model. How might the susceptible population be replenished?

        Over short time horizons those who have been previously infected may lose their immunity allowing the potential for reinfection. For COVID-19, estimates for immunity duration are just a few months. Over longer time horizons an additional source of susceptible individuals may be via population demography with births adding to the population. Both of these processes can be incorporated to extend the SIR model. 
        Here, we will focus on adding waning immunity. Such models are often termed SIRS to note that previously infected individuals may be susceptible to reinfection in the future.

        Our amended model can be written as:
        $$\frac{dS}{dt} = \overset{transmission}{\overbrace{-βSI}} + \overset{waning\ immunity}{\overbrace{αR}}$$

        $$\frac{dI}{dt} = \overset{transmission}{\overbrace{βSI}} - \overset{recovery}{\overbrace{γI}}$$

        $$\frac{dR}{dt} = \overset{recovery}{\overbrace{γI}} - \overset{waning\ immunity}{\overbrace{αR}}.$$

        There are only small changes from the SIR model here! Note, as before, γ represents the rate for an infectious individual to stop being infectious (such that the average recovery time is 1/γ). Now, we use α to represent the rate of waning immunity, such that the average time to move from R back to S is 1/α. All together this suggests the average time that an individual who has just been infected before being susceptible again can be computed as: 1/γ + 1/α.
        """
    )
    return


@app.cell
def _(beta2, gamma2, inversealpha, odeint, t2, u0):
    # The SIRS model differential equations.
    def SIRSmodel(u, t, beta, gamma, alpha):
        S, I, R = u
        dSdt = -beta * S * I + alpha * R
        dIdt = beta * S * I - gamma * I
        dRdt = gamma * I - alpha * R
        return dSdt, dIdt, dRdt

    ret2 = odeint(SIRSmodel, u0, t2, args=(beta2, gamma2, 1/(inversealpha.value)))
    S2, I2, R2 = ret2.T
    return I2, R2, S2, SIRSmodel, ret2


@app.cell(hide_code=True)
def _(mo, np):
    beta2 = 0.4
    gamma2 = 0.1
    t2 =  np.linspace(0, 1000,1000)
    inversealpha = mo.ui.slider(1, 1000, value=365, step = 1, label='Average immunity duration 1/α  (days): ')
    sim_length = mo.ui.dropdown(options={"200 days":200, "3 years":1040},
                            value="200 days", # initial value
                            label="Choose length to time to view simulation output for")
    mo.md(
        f"""
        **Defining epidemiological parameters for the SIRS model.**

        Transmission rate, β = 0.4 /day    (note, this is fixed)

        Recovery rate, γ = 0.1/day    (note, this is fixed)

        {inversealpha}

        {sim_length}
        """
    )
    return beta2, gamma2, inversealpha, sim_length, t2


@app.cell
def _(I2, R2, S2, plt, sim_length, t2):
    plt.figure(figsize=(8,4))
    plt.plot(t2,S2,label='S')
    plt.plot(t2,I2,label='I')
    plt.plot(t2,R2,label='R')
    plt.xlabel('Time (days)')
    plt.ylabel('Population fraction')
    plt.xlim([-10,sim_length.value+5])
    plt.legend()
    plt.title("Simulation of epidemiological population dynamics with the SIRS model")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        ### Question:
        What are the differences in dynamics between the SIR and the SIRS model? Hint: go back to the SIR model and match the parameters used here in the SIRS model output.
        """
    ).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        # 4) The stochastic SIR model
        The differential equation framework utilized by the SIR and SIRS models above provides a 'mean-field' description of the systems dynamics -- that is, they provide dynamical descriptions that assume that populations have average behavior. This works well when population sizes are large, but there may be deviations when population sizes are smaller due to randomness a.k.a stochasticity. To explore this we will revisit the SIR model, but use a different framework to consider how the dynamics change.
        In this new framework we will consider individual 'events' i.e., we will track each individual infection and each recovery in the population, and rather than population fractions we will consider a population of a particular size. 

        **Q:** what happens during an infection event? **A:** The number of susceptibles decreases by 1, and the number of infectious individuals increases by 1.   So, what happens during a recovery event?

        In the simulation, we will proceed by randomly choosing one event at a time - but not with 50:50 odds. When there are many susceptible individuals we might expect that transmission events are more likely than recovery events; and when there are many infectious individuals the chance of a recovery event occuring will increase. The chance of the next event being a transmission or recovery depends on the state of the system, the current number of susceptible, infectious and recovered individuals.

        Note that the time between events will also differ depending on the state of the system, when few people are infected we may expect a longer time for the next event vs. when many people are infected.

        The simulation technique we use here is termed the Gillespie algorithm (aka stochastic simulation algorithm): https://en.wikipedia.org/wiki/Gillespie_algorithm

        Our implementation of the Gillespie algorithm for the SIR model is in the next code block:
        """
         )
    return


@app.cell
def _(Pop, np):
    import random #allows us to generate random numbers --- random.random() generates a random number between 0 and 1.

    def stochasticSIR(InitS,InitI,InitR,beta,gamma):
        #initialize some lists to store and append values as we step from one event to the next
        time = [] #list to save times
        Slist = [] #list to save S
        Ilist = [] #list to save I
        Rlist = [] #list to save R
        time.append(0) #at intial time 0
        Slist.append(InitS) #at intial time 0 S will have SInit individuals
        Ilist.append(InitI) #at intial time 0 I will have IInit individuals
        Rlist.append(InitR) #at intial time 0 R will have RInit individuals
        Itot=InitI
        #start the reaction loop
        while Itot>0: #while more than one infectious individual in the population
            r = np.array([beta*Slist[-1]*Ilist[-1]/Pop.value, gamma*Ilist[-1]]) #rates of transmission and recovery in the population
            if sum(r<0)>0: # a check to make sure rates are not negative, and if yes to exit.
                break
            sumr = sum(r)    #total rates at this point in time in the population
            r_prob = np.cumsum(r)/sumr;
            tau = 1/sumr *np.log(1/random.random()); #time to next event
            time.append(time[-1]+tau); #access the last element in the current list of time, and add tau
            choose = random.random() #random value between 0 and 1
            index = np.where(choose>r_prob)[0] #is the random value greater than the proportional rates?
            if len(index)==1: 
                Itot = Itot-1 #recovery event
                Slist.append(Slist[-1]) 
                Ilist.append(Ilist[-1] - 1) 
                Rlist.append(Rlist[-1] + 1)
            else:
                Itot = Itot+1 #transmission event
                Slist.append(Slist[-1] - 1)
                Ilist.append(Ilist[-1] + 1)
                Rlist.append(Rlist[-1])

        return time, Slist, Ilist, Rlist
    return random, stochasticSIR


@app.cell(hide_code=True)
def _(mo):
    Pop = mo.ui.dropdown(options={"10,000":10000, "500":500, "100":100},value="10,000", 
                            label="choose population size")#population size
    buttonsirstoch = mo.ui.run_button()


    mo.md(
        f"""
        **Epidemiological parameters for the stochastic SIR model.**

        Transmission rate, β = 0.4 /day 

        Recovery rate, γ = 0.1/day

        {Pop}

        When you are ready press this button to run the simulation:
        {buttonsirstoch}
        """
    )
    return Pop, buttonsirstoch


@app.cell
def _(Pop, buttonsirstoch, mo, stochasticSIR):
    mo.stop(not buttonsirstoch.value) #don't run this code until the button above is pressed.

    #inital conditions for simulation
    InitI = 1 # initial infectious individuals in the population
    InitR = 0 # initial recovered individuals in the population
    InitS = Pop.value - InitI - InitR # initial susceptibile individuals in the population

    beta3 = 0.4 #transmission rate for stochastic SIR simulations
    gamma3 = 0.1 #recovery rate for stochastic SIR simulations

    ret3 = stochasticSIR(InitS,InitI,InitR,beta3,gamma3)
    time, Slist, Ilist, Rlist = ret3

    mo.md(
        f"""
        Now we will initialze the simulation with one infected individual in the population and run the dynamics, which will show in the plot below.

        **Note:** try pressing the 'click to run' button above again. Do you get the same result?
        """
    )
    return (
        Ilist,
        InitI,
        InitR,
        InitS,
        Rlist,
        Slist,
        beta3,
        gamma3,
        ret3,
        time,
    )


@app.cell
def _(Ilist, Rlist, Slist, plt, time):
    plt.figure(figsize=(8,4))
    plt.plot(time,Slist,label="S")
    plt.plot(time,Ilist,label="I")
    plt.plot(time,Rlist,label="R")
    plt.legend()
    plt.xlabel("Time (days)")
    plt.ylabel("People")
    plt.title("Stochastic SIR model simulation")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        ### Question
        Do the dynamics match up to the SIR model in section 1? Are they similar in their qualitative behavior, and the timing of epidemic peaks? How does this change with population size? 
        """
    ).callout("success")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
            ### Did you spot this?
            As we initialize this simulation with just one infected individual, there is a chance that the epidemic will fizzle out (quickly reduce number of infectious individuals to 0) rather than lead to a full blown outbreak due to stochasticity. In particular, the outbreak probability for the SIR model can be calculated as:
            $$P = 1 - \left(\frac{1}{\mathcal{R}_0}\right)^m$$
            where m is the number of initially infected individuals. The time horizon in simulations that fizzle out may be (very) short. Even though the epidemiological context is the same, each of our model simulations here shows a different trajectory! This is *just one* type of uncertainty that can make it difficult to anticipate and forecast how an epidemic may change. 

            **Q: What is the outbreak probability for the stochastic model we simulated above?**
            """
        ).callout("info")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        #5) Viewing the distribution of outbreak sizes in the stochastic SIR model

        Note: Press the 'click to run' button below to generate a plot of outbreak sizes calculated by running the stochastic SIR model above 100 times. Do be aware that this can take a while to compute when the population size is large -- many events can occur when the population is large. Consider, first going back up and changing the population size.

        Compare with your peers! Do the distributions match?
        """
        )
    return


@app.cell(hide_code=True)
def _(InitI, InitR, InitS, beta3, button, gamma3, mo, stochasticSIR):
    mo.stop(not button.value) # don't run the below code until the button is pressed.
    listofR = []
    #run the stochastic SIR model 100 times
    for ii in range(100): #repeat this code while ii takes values from 0 to 99.
        ret4 = stochasticSIR(InitS,InitI,InitR,beta3,gamma3)
        listofR.append(ret4[3][-1]) #add the number of recovered individuals at the end of each simulation to the list listofR.
    return ii, listofR, ret4


@app.cell(hide_code=True)
def _(mo):
    button = mo.ui.run_button()
    button
    return (button,)


@app.cell
def _(listofR, np, plt):
    plt.figure(figsize=(8,4))
    plt.hist(np.array(listofR),bins=100)
    plt.xlabel("Outbreak Size")
    plt.ylabel("Number of outbreaks")
    plt.title("Outbreak size distributions from 100 simulations")
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        # You made it!
        Congratulations you made it to the end of the notebook tutorial! 👏👏👏
        If there is still time, think about what these simple models say...and what they are missing!
        If you are feeling adventurous, consider going back through the code and consider changing values e.g. in cell 17 try running the SIRS model with different values for the transmission rate and recovery rate.
        """
    ).callout("success")
    return


if __name__ == "__main__":
    app.run()
