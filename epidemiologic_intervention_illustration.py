"""
Kucharavy Andrei 2020 all rights reserved - MIT license.

Please note that this simulation engine does not make an epidemiologist out of you.
Much more work goes into simulating and calculating things properly

You will need a base python installation to run this.

We use the following assumptions about the disease:
 - 5 days, non-infective incubation period
 - 48 hours pre-symptomatic infectious period
 - 14 days symptomatic infections period
 - after symptoms subside, the population acquires immunity
All of the above can be adjusted using the `state_duration` list

We use the following assumptions about the transmission:
 - any points within a close circle ("sure contaminations") will be contaminated as soon as an
 individual becomes infectious. It can be close coworkers/household members/hobby activity contacts
 - at most 1 individual in a larger circle ("possible contaminations") can be contaminated per
 day. The chance of contamination is proportional to the distance to the contagious individual.
 There is a chance that on a given day no contamination occurs that is set as a parameter.
All of the above can be adjusted using the `limit_distances` list

We use the following assumptions about the measures:
 - The quarantine is not perfect, but reduces the radius of "sure contaminations" and "possible
 contaminations" circles as well as raises the probability no contamination occurs on a given day in
 "possible contaminations" circle (eg transmissions within the household, during store visits, ...)
 - The degree of efficiency of quarantine is controlled by the `quarantine power` parameter
 - The adherence to all the measures (quarantine, self-isolation, contact tracing & isolation...) is
 total
 - Contact tracing can vary between two states: "from memory only", where only the "sure
 contaminations" are isolated and "total tracking", where in addition all the "possible
 contaminations" are isolated as well.
 - The degree of efficiency of the contact tracing power is controlled by the
 `contact_tracing_power` parameter.

R_t is calculated as the average number of individuals infected by a single infectious
individual across all the individuals that are still infections

"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import spatial
import random
from matplotlib import colors as mcol
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


# repeatability of simulations
random.seed(1981)
np.random.seed(1981)


def run_simulation():

    points_x = np.random.uniform(size=(individuals,))
    points_y = np.random.uniform(size=(individuals,))

    combines_positions = np.vstack((points_x, points_y)).T
    dist_matrix = spatial.distance.cdist(combines_positions, combines_positions)

    naive = 0.
    incubation = 0.75
    pre_sympt_inf = 0.5
    sympt_inf = 0.375
    immune = 1.

    infection_status = np.random.choice([immune, naive],
                                        size=(individuals,),
                                        p=[fraction_initially_immune, 1 - fraction_initially_immune])

    immune_from_start = np.sum(infection_status)

    isolation_status = np.zeros_like(infection_status).astype(np.bool)

    duration = - np.ones_like(infection_status)

    fig, ax = plt.subplots()

    ax1 = plt.subplot(212)
    scatter = ax1.scatter(points_x, points_y, c=infection_status, cmap='Dark2')

    ax2 = plt.subplot(221)
    ax2.set_title('R_t')
    R_t_plot = ax2.plot([], [], '-ok')
    ax2.set_xlabel('days since first infection')
    ax2.set_ylabel('R_t')
    # ax2.set_ylim(bottom=0)

    ax3 = plt.subplot(222)
    ax3.set_title('population_state')
    active_cmap = plt.get_cmap('Dark2')
    naive_plt = ax3.plot([], [], '-o', c=active_cmap(naive))
    incub_plt = ax3.plot([], [], '-o', c=active_cmap(incubation))
    pre_symp_plt = ax3.plot([], [], '-o', c=active_cmap(pre_sympt_inf))
    sympt_plt = ax3.plot([], [], '-o', c=active_cmap(sympt_inf))
    immune_plt = ax3.plot([], [], '-o', c=active_cmap(immune))
    ax3.set_yscale('log')
    ax3.set_xlabel('days since first infection')
    ax3.set_ylabel('log individuals in the state x')

    fig.set_size_inches(19, 10.5)

    print('fig size: {0} DPI, size in inches {1}'.format(
        fig.get_dpi(), fig.get_size_inches()))

    sub_cmap = plt.get_cmap('Set1')
    quarantine_color = sub_cmap(0)
    self_isolation_color = sub_cmap(0.45)
    contact_tracing_color = sub_cmap(0.23)


    patient_zero = None
    R_t_spread_to = np.zeros_like(infection_status)
    R_t_tracker = []
    population_state_tracker = []
    quarantine = False
    self_isolation = False
    contact_tracing = False
    simulation_day = 0
    total_cases = 0
    max_symptomatic = 0

    def update(i):
        global limit_distances
        nonlocal patient_zero, quarantine, simulation_day, isolation_status, \
            self_isolation, contact_tracing, total_cases, max_symptomatic

        quarantine_limit_distances = (limit_distances[0] / quarantine_power,
                                        limit_distances[1] / quarantine_power,
                                        limit_distances[2] * quarantine_power)

        simulation_day += 1
        i = simulation_day

        ax_2_q_traced = True
        ax_3_q_traced = True

        ax_2_qe_traced = True
        ax_3_qe_traced = True

        ax_2_si_traced = True
        ax_3_si_traced = True

        ax_2_ct_traced = True
        ax_3_ct_traced = True


        try:
            if i == quarantine_day:
                # print('quarantine_today')
                quarantine = True
                isolation_status[:] = True
                ax_2_q_traced = False
                ax_3_q_traced = False

            if i == quarantine_end_day:
                quarantine = True
                isolation_status[:] = False
                ax_2_qe_traced = False
                ax_3_qe_traced = False

            if i == self_isolation_day:
                self_isolation = True
                ax_2_si_traced = False
                ax_3_si_traced = False

            if i == contact_tracing_day:
                contact_tracing = True
                ax_2_ct_traced = False
                ax_3_ct_traced = False

            if i == 0:
                active_spreaders = []
                pass

            elif i == 1:
                patient_zero = np.random.choice(np.arange(0, individuals)[infection_status == naive])
                # print('patient zero:', patient_zero)
                infection_status[patient_zero] = incubation
                scatter.set_array(infection_status)
                duration[patient_zero] = 0
                active_spreaders = [patient_zero]

            else:
                counter_on = np.arange(0, individuals)[np.logical_and(
                    infection_status != naive, infection_status != immune)]
                duration[counter_on] += 1
                infection_status[duration >= 0] = incubation
                infection_status[duration >= state_duration[1]] = pre_sympt_inf
                infection_status[duration >= state_duration[2]] = sympt_inf
                infection_status[duration >= state_duration[3]] = immune
                duration[duration >= state_duration[3]] = -1

                if contact_tracing:
                    roots_of_isolation = np.arange(individuals)[
                        np.logical_and(infection_status == sympt_inf,
                                       np.logical_not(isolation_status))]

                    # print(roots_of_isolation)

                    for root in roots_of_isolation:
                        # print('processing root', root)

                        reach_array = dist_matrix[root, :]

                        potential_infections_pad = reach_array < limit_distances[1]

                        certain_infections_pad = reach_array < limit_distances[0]

                        # isolation_status[certain_infections_pad] = True

                        potential_infections_pad[certain_infections_pad] = False

                        # print('mapping complete',
                        #       np.sum(potential_infections_pad),
                        #       np.sum(certain_infections_pad))

                        if np.any(certain_infections_pad):
                           to_isolate_surely = np.arange(0, individuals)[certain_infections_pad]
                           # print('sure isolation: ', to_isolate_surely)
                           isolation_status[to_isolate_surely] = True

                           # isolation_status[potential_infections_pad] = True

                        if np.any(potential_infections_pad):
                            corrected_reach_array = reach_array[potential_infections_pad]
                            corrected_reach_array = corrected_reach_array / np.sum(corrected_reach_array)
                            reach = int(np.ceil(np.sum(
                                potential_infections_pad)*contact_tracing_power))
                            # print('reach', reach, np.sum(np.sum(potential_infections_pad)))
                            # if reach == np.sum(np.sum(potential_infections_pad)):
                            #     isolation_status[potential_infections_pad] = True

                            if reach > 0:
                                to_isolate = np.random.choice(
                                    np.arange(0, individuals)[potential_infections_pad],
                                    (reach,), replace=False,
                                    p=corrected_reach_array)
                                # print(set(np.arange(0, individuals)[
                                #               potential_infections_pad].tolist()) -
                                #       set(to_isolate.tolist()))
                                isolation_status[to_isolate] = True

                if self_isolation:
                    isolation_status[infection_status == sympt_inf] = True
                    # print(isolation_status)


                # print('infection_status', infection_status[np.logical_and(
                #     infection_status != naive, infection_status != immune)])
                # print('duration', duration[duration != -1])

                active_spreaders = np.arange(0, individuals)[
                        np.logical_or(infection_status == pre_sympt_inf,
                                      infection_status == sympt_inf)]

                # print('active spreaders', active_spreaders)

                for spreader in active_spreaders:
                    reach_array = dist_matrix[spreader, :]

                    if isolation_status[spreader]:
                        _limit_distances = quarantine_limit_distances
                    else:
                        _limit_distances = limit_distances

                    potential_infections_pad = np.logical_and(reach_array < _limit_distances[1],
                                                              infection_status == naive)

                    certain_infections_pad = np.logical_and(reach_array < _limit_distances[0],
                                                              infection_status == naive)

                    potential_infections_pad[certain_infections_pad] = False

                    if np.any(certain_infections_pad):
                        infected_nodes = np.arange(0, individuals)[certain_infections_pad]
                        # print('certainly infected:', infected_nodes)
                        infection_status[infected_nodes] = incubation

                        for infected in infected_nodes:
                            R_t_spread_to[spreader] += 1
                            # print('plotting spread to', infected)
                            ax1.plot([points_x[spreader], points_x[infected]],
                                     [points_y[spreader], points_y[infected]],
                                     '-r', zorder=0)

                    if np.any(potential_infections_pad) and np.random.uniform() > \
                            _limit_distances[2]:
                        corrected_reach_array = reach_array[potential_infections_pad]
                        corrected_reach_array = corrected_reach_array / np.sum(corrected_reach_array)
                        infected = np.random.choice(
                            np.arange(0, individuals)[potential_infections_pad],
                            p=corrected_reach_array)
                        infection_status[infected] = incubation
                        R_t_spread_to[spreader] += 1
                        ax1.plot([points_x[spreader], points_x[infected]],
                                 [points_y[spreader], points_y[infected]],
                                 '-r', zorder=0)

                # print('trackeur R_t', R_t_spread_to[R_t_spread_to > 0])

                scatter.set_array(infection_status)
                ring_color = [active_cmap(val) for val in infection_status]

                for _i in np.arange(individuals)[isolation_status]:
                    ring_color[_i] = mcol.to_rgb('k')

                scatter.set_edgecolor(ring_color)

            try:
                R_t = np.mean(R_t_spread_to[active_spreaders])
                R_t_tracker.append(R_t)

                if i > 12:
                    ax2.plot([i-1, i], [R_t_tracker[-2], R_t_tracker[-1]], '-ok')

                if not ax_2_q_traced:
                    ax2.axvline(quarantine_day, color=quarantine_color)
                    ax_2_q_traced = True
                if not ax_2_qe_traced:
                    _, y_top = ax2.get_ylim()
                    q_rec = Rectangle((quarantine_day, 0),
                                      quarantine_end_day-quarantine_day, y_top)
                    pc = PatchCollection([q_rec],  alpha=0.3,
                                         facecolors=quarantine_color,
                                         edgecolors=quarantine_color)
                    ax2.add_collection(pc)
                    ax2.axvline(quarantine_end_day, color=quarantine_color)
                    ax_2_qe_traced = True
                if not ax_2_si_traced:
                    ax2.axvline(self_isolation_day, color=self_isolation_color)
                    ax_2_si_traced = True
                if not ax_2_ct_traced:
                    ax2.axvline(contact_tracing_day, color=contact_tracing_color)
                    ax_2_ct_traced = True


                population_state_tracker.append([np.sum(infection_status == naive),
                                                np.sum(infection_status == incubation),
                                                np.sum(infection_status == pre_sympt_inf),
                                                np.sum(infection_status == sympt_inf),
                                                np.sum(infection_status == immune)])

                drawable_pop_state = np.array(population_state_tracker)
                # print(drawable_pop_state.shape)
                # print(drawable_pop_state[:, 0].shape)

                if i > 1:
                    ax3.plot([i-1, i],
                             drawable_pop_state[-2:, 0],
                             '-o', c=active_cmap(naive), label='naive')
                    ax3.plot([i-1, i],
                             drawable_pop_state[-2:, 1],
                             '-o', c=active_cmap(incubation), label='incub')
                    ax3.plot([i-1, i],
                             drawable_pop_state[-2:, 2],
                             '-o', c=active_cmap(pre_sympt_inf), label='pre-sympt')
                    ax3.plot([i-1, i],
                             drawable_pop_state[-2:, 3],
                             '-o', c=active_cmap(sympt_inf), label='symptomatic')
                    ax3.plot([i-1, i],
                             drawable_pop_state[-2:, 4],
                             '-o', c=active_cmap(immune), label='immune')
                if i == 2:
                    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

                if not ax_3_q_traced:
                    ax3.axvline(quarantine_day, color=quarantine_color)
                    ax_3_q_traced = True
                if not ax_3_qe_traced:
                    q_rec = Rectangle((quarantine_day, 0),
                                      quarantine_end_day-quarantine_day, individuals)
                    pc = PatchCollection([q_rec],  alpha=0.3,
                                         facecolors=quarantine_color,
                                         edgecolors=quarantine_color)
                    ax3.add_collection(pc)
                    ax3.axvline(quarantine_end_day, color=quarantine_color)
                    ax_3_qe_traced = True
                if not ax_3_si_traced:
                    ax3.axvline(self_isolation_day, color=self_isolation_color)
                    ax_3_si_traced = True
                if not ax_3_ct_traced:
                    ax3.axvline(contact_tracing_day, color=contact_tracing_color)
                    ax_3_ct_traced = True

            except Exception as exc:
                print(exc)

            description_string = ''
            if quarantine:
                description_string += ', quarantine'
            if self_isolation:
                description_string += ', self-isolation'
            if contact_tracing:
                description_string += ', contact-tracing @ %.2f%%' % (contact_tracing_power*100)

            label = 'day %d: R_t: %.2f, %s\n' \
                        'naive: %d, incub: %d, pre-sympt :%d ' \
                        'sympt: %d, immune: %d' % (i, R_t, description_string,
                                                   *population_state_tracker[-1])

            ax1.set_xlabel(label)
            print(label)

        except Exception as exc:
            print(exc)

        if population_state_tracker[-1][1] + population_state_tracker[-1][2] \
                    + population_state_tracker[-1][3] == 0:

            total_cases = population_state_tracker[-1][-1] - immune_from_start
            max_symptomatic = np.max(drawable_pop_state[:, 3])

            print('Epidemy stopped on day %d with %d cases' % (i, total_cases))
            raise Exception('Epidemy stopped on day %d with %d cases' % (i, total_cases))

        return scatter, ax, R_t_plot, naive_plt, incub_plt, pre_symp_plt, sympt_plt, immune_plt

    anim = FuncAnimation(fig, update, 1000, interval=500)

    name_payload = []
    if quarantine_day > 0:
        name_payload.append('quarantine_power_%d_starting_on_%d' %
                            (quarantine_power, quarantine_day))
        if quarantine_end_day > 0:
            name_payload.append('ending_on_%d' % quarantine_end_day)

    if self_isolation_day > 0:
        name_payload.append('self_isolation_on_%d' % self_isolation_day)

    if contact_tracing_day > 0:
        name_payload.append('contact_tracing_on_%d_with_power%.2f' %
                            (contact_tracing_day, contact_tracing_power))

    if len(name_payload) == 0:
        name_payload.append('base_scenario')

    name_payload.append('contagion_parameters_%.2f_%.2f_%.2f' % limit_distances)

    name_payload.append('stopped_on_%d_with_%d_cases' % (simulation_day, total_cases))
    name_payload.append('%d_max_symptimatic_cases' % max_symptomatic)

    name_struct = '-'.join(name_payload)

    name_struct += '.gif'

    if save:
        anim.save(name_struct, dpi=80, writer='imagemagick')
    else:
        plt.show()


if __name__ == '__main__':
    # Assumption about stages of disease duration.
    # All durations is from the contamination day (day 0). -1 means infinite duration.
    # [naive, incubation, pre-symptomatic infectious, symptomatic infectious, immune]
    state_duration = [-1, 5, 7, 21, -1]

    # Total individuals in the simulation
    individuals = 1000

    # Euclidian distances on a plane that a disease will be able to jump to contaminate
    # (sure contamination, possible contamination, chance not to contaminate anyone in a
    # "possible circle" on a given day)
    limit_distances = (0.05, 0.1, 0.1)

    # Fraction of individuals initially immune. Due to the pecularities of how colors work in
    # matplotlib's scatterplot, at least one individual needs to be immune at the beginning.
    # An alternative interpretation is the fraction of silent, asymptomatic and non-contagious
    # infections
    fraction_initially_immune = 0.1

    # For all of the parameters below, the number is the day of the simulation at which the
    # measaure is taken. -1 means the measure is never taken
    quarantine_day = 20  # day at which a universal quarantine is put in place
    quarantine_end_day = 40  # day at which the universal quarantine ends
    self_isolation_day = 30  # day at which everyone with symptoms start to self-quarantine
    contact_tracing_day = 35  # day at which the contact tracing and contact quarantining starts

    # Parameters below decide how well the measures put in place are implemented.

    # Quarantine power is the divider of the parameters of transmission (decreases the closest
    # "sure contamination" and "possible contamination" circles radiuses, by dividing the base
    # amount by the `quarantine_power` value )
    quarantine_power = 4

    # Contact tracing power is how well the contacts are traced.
    # 0. means only the "sure contamination" contacts - likely those coming from work and home,
    # that the person can inform as soon as they become symptomatic - are tracked and put into
    # quarantine.
    # 1. means all the contacts in "sure contaminations" and "possible contaminations" circles
    # are tracked down and put into quarantine (eg by an app)
    # Values in between the two are the fraction of "possible contaminations" that are tracked
    contact_tracing_power = 1.

    # Flag as to whether save the results or just show in a pyplot window.
    save = True

    # command that actually runs the simulation
    run_simulation()

    # reference durations:

    # nothing: ~ 106 days
    # quarantine to 4: 63
    # quarantine to 2 on 20, for 211 days
    # quarantine at 2 on 20, for 6 weeks: 131 days
    # self-isolation on 20, power 4: 106 days
    # self-isolation and total tracking on 20, power 4: 77 days
    # self-isolation and partial tracking on 20, power 4: 122 days
