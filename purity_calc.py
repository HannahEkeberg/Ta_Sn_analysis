from get_variables import *
import curie as ci

class purity_calc:

    def __init__(self, foil, isotopes): #'Ta01'
        stack_numbs = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12', '13', '14']
        foils = ['Ta' + number for number in stack_numbs]
        self.foil = foil
        self.foil_idx = foils.index(foil)
        self.isotopes = isotopes
        # self.isotope = isotope
        # self.A0_177W = eob_activity_manually(element='Ta', isotope='177W', independent=True, stack=None)[0]
        # self.A0_178W = eob_activity_manually(element='Ta', isotope='178W', independent=None, stack=None)[0]
        # self.A0_179W = eob_activity_manually(element='Ta', isotope='179W', independent=False, stack=None)[0]
        # self.A0_181W = eob_activity_manually(element='Ta', isotope='181W', independent=True, stack=None)[0]
        # self.A0_176W = eob_activity_manually(element='Ta', isotope='176W', independent=True, stack=None)[0]

    # def get_purity(self):
    #     cooling_times  = np.linspace(0,3600*20, 1000) # hours
    #     # decay_constant_177W = ci.Isotope('177Wg').decay_const()




    #     A0_177W = self.A0_177W[self.foil_idx]
    #     A_177W = self.getA(A0_177W, decay_constant_177W, cooling_times)

    #     # A_178W = self.getA(self.A0_178W[self.foil_idx], ci.Isotope('178Wg').decay_const(), cooling_times)    
    #     # A_179W = self.getA(self.A0_179W[self.foil_idx], ci.Isotope('179Wg').decay_const(), cooling_times)    
    #     # A_181W = self.getA(self.A0_181W[self.foil_idx], ci.Isotope('181Wg').decay_const(), cooling_times)    
    #     # A_176W = self.getA(self.A0_176W[self.foil_idx], ci.Isotope('176Wg').decay_const(), cooling_times)   


    #     # total_activity = A_177W +  A_178W + A_179W + A_176W# + A_181W


    #     if A_177W[0]>1e-4:
    #         yield_curve = A_177W/total_activity*100
    #         idx = np.argmax(yield_curve)
    #         max_time = cooling_times[idx] / 3600
    #         max_yield = yield_curve[idx]
    #         plt.axvline(max_time, linestyle='--', linewidth=0.5,
    #         label=f'Max yield: {max_time:.1f} h ( {max_yield:.1f} %)')
    #         print(f'Max yield: {max_time:.1f} h ( {max_yield:.1f} %)')
    #         plt.plot(cooling_times/3600, yield_curve, label='purity 177W', color='#741b47')
    #     if A_178W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_178W/total_activity*100, label='purity 178W', color='dodgerblue', linewidth=1.2, linestyle='-.')
    #     if A_179W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_179W/total_activity*100, label='purity 179W', color='#cf9b65', linewidth=1.2, linestyle='-.')
    #     # plt.plot(cooling_times/3600, A_181W/total_activity*100, label='purity 181W', color='blue')
    #     if A_176W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_176W/total_activity*100, label='purity 176W', color='#ac6186', linewidth=1.2, linestyle='-.')
    #     plt.legend()
    #     plt.xlabel('Hours since end of beam')
    #     plt.ylabel('Radioisotopic purity (%)')
    #     plt.title(self.foil)
    #     plt.show()
    #     # plt.savefig('generatedfiles/yields/' + self.foil + '.pdf')
    #     # plt.close()

    def purity(self, isotope):
        cooling_times  = np.linspace(0,3600*20, 1000) # hours
        yield_curve = self.yield_curve(isotope, cooling_times)
        return cooling_times, yield_curve


    def print_yield(self, yield_curve, cooling_times):
        idx = np.argmax(yield_curve)
        max_time = cooling_times[idx] / 3600
        max_yield = yield_curve[idx]
        print(f'Max yield: {max_time:.1f} h ( {max_yield:.1f} %)')

    def yield_curve(self, isotope, t):
        A0 = self.A0(isotope)[self.foil_idx]
        A = self.getA(A0, isotope, t)
        total_activity = self.total_activity(t)
        return A/total_activity * 100

    def total_activity(self, t):
        eob_activity = []
        activity = []
        for i in self.isotopes:
            print(i)
            A0 = self.A0(i)[self.foil_idx]
            if A0 > 1e-4:
                eob_activity.append(A0)
                activity.append(self.getA(A0, i, t))
        print(len(activity))
        total = np.sum(activity, axis=0) #np.sum(eob_activity) if len(eob_activity)>0 else 1
        return total

    def A0(self, isotope):
        try:
            return eob_activity_manually(element='Ta', isotope=isotope, independent=True, stack=None)[0]
        except: 
            return eob_activity_manually(element='Ta', isotope=isotope, independent=None, stack=None)[0]
        # return eob_activity_manually(element='Ta', isotope=isotope, independent=True, stack=None)[0]
        # return eob_activity_manually(element='Ta', isotope=isotope, independent=True, stack=None)[0]
    

    def getA(self, A0, isotope, t):
        decay_constant = ci.Isotope(isotope).decay_const()
        return A0 * np.exp(-decay_constant * t)


def plot():
    pc = purity_calc('Ta01', ['176', '178', '177'])
    W178 = pc.purity('178W') 
    plt.plot(W178[0]/3600, W178[1], label='purity 177W', color='#741b47')
    plt.show()
    #     if A_178W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_178W/total_activity*100, label='purity 178W', color='dodgerblue', linewidth=1.2, linestyle='-.')
    #     if A_179W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_179W/total_activity*100, label='purity 179W', color='#cf9b65', linewidth=1.2, linestyle='-.')
    #     # plt.plot(cooling_times/3600, A_181W/total_activity*100, label='purity 181W', color='blue')
    #     if A_176W[0]>1e-4:
    #         plt.plot(cooling_times/3600, A_176W/total_activity*100, label='purity 176W', color='#ac6186', linewidth=1.2, linestyle='-.')
    #     plt.legend()
    #     plt.xlabel('Hours since end of beam')
    #     plt.ylabel('Radioisotopic purity (%)')
    #     plt.title(self.foil)
    #     plt.show()
    #     # plt.savefig('generatedfiles/yields/' + self.foil + '.pdf')
    #     # plt.close()


plot()

# t1,y1 = purity_calc('Ta01', ['176W', '178W', '179W']).purity('177W')
# t2,y2 = purity_calc('Ta01', ['178W', '177W',  '179W']).purity('176W')
# t3,y3 = purity_calc('Ta01', ['176W', '177W', '179W']).purity('178W')
# # t4,y4 = purity_calc('Ta01', ['176W', '177W', '178W']).purity('179W')
# plt.plot(t1/3600,y1, label='1')
# plt.plot(t2/3600,y2, label='2')
# plt.plot(t3/3600,y3, label='3')
# # plt.plot(t4/3600,y4, label='4')
# plt.legend()
# plt.show()
# purity_calc('Ta02').get_purity()
# purity_calc('Ta03').get_purity()
# purity_calc('Ta04').get_purity()
# purity_calc('Ta05').get_purity()
# purity_calc('Ta06').get_purity()
# purity_calc('Ta07').get_purity()



# purity_calc('Ta06').get_purity()
# purity_calc('Ta07').get_purity()
# purity_calc('Ta08').get_purity()
# purity_calc('Ta09').get_purity()
# purity_calc('Ta10').get_purity()
# purity_calc('Ta11').get_purity()
# purity_calc('Ta12').get_purity()
# purity_calc('Ta13').get_purity()
# purity_calc('Ta14').get_purity()