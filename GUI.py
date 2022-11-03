import tkinter
import tkinter.messagebox
import customtkinter
import gmpy2
import primefinder
import random
import TestV2
from sympy.ntheory import factorint
from PIL import Image, ImageTk
customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

#global k
#k = 0

class App(customtkinter.CTk):

    WIDTH = 780
    HEIGHT = 520
   # Ind = 0
    filter = True
    canonical = True
    primesize = 0
    i = 0
    prime =0
    A = 0
    B=0
    order = "0"
    userset = True
    method = 2
    def __init__(self):
        super().__init__()
        #self.resizable(False,False)
        self.title("Elliptic Curve Arithmetic")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # call .on_closing() when app gets closed

        # ============ create two frames ============

        # configure grid layout (2x1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self,
                                                 width=180,
                                                 corner_radius=0)

        self.frame_left.grid(row=0, column=0, sticky="nswe")
        self.frame_left.grid_propagate(False)
        self.frame_right = customtkinter.CTkFrame(master=self)
        self.frame_right.grid(row=0, column=1, sticky="nswe", padx=20, pady=20)

        # ============ frame_left ============

        # configure grid layout (1x11)
        self.frame_left.grid_rowconfigure(0, minsize=10)   # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(5, weight=1)  # empty row as spacing
        self.frame_left.grid_rowconfigure(8, minsize=20)    # empty row with minsize as spacing
        self.frame_left.grid_rowconfigure(11, minsize=10)  # empty row with minsize as spacing
        self.img = ImageTk.PhotoImage(Image.open("v2-removebg-preview.png"))
        self.label_1 = customtkinter.CTkLabel(master=self.frame_left,text="Tools",
                                              text_font=("Roboto Medium", -16))  # font name and size in px
        self.label_1.grid(row=1, column=0, pady=10, padx=10)


        self.button_2 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Point Counting"
                                                )
        self.button_2.grid(row=3, column=0, pady=10, padx=20)
        self.button_8 = customtkinter.CTkButton(master=self.frame_left,
                                                text="ECDLP",
                                                command=self.button_event1)
        self.button_8.grid(row=4, column=0, pady=10, padx=20)
        self.button_9 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Pairings",
                                                command=self.button_event1)
        self.button_9.grid(row=5, column=0, pady=10, padx=20)
        self.button_10 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Elliptic Curve Factorisation",
                                                command=self.button_event1)
        self.button_10.grid(row=6, column=0, pady=10, padx=20)
        self.button_11 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Basic Arithmetic",
                                                command=self.button_event1)
        self.button_11.grid(row=7, column=0, pady=10, padx=20)
        self.button_12 = customtkinter.CTkButton(master=self.frame_left,
                                                text="Help",
                                                command=self.button_event1)
        self.button_12.grid(row=8, column=0, pady=10, padx=20)

        self.button_3 = customtkinter.CTkButton(master=self.frame_right,
                                                text="Generate random curve",
                                                 width=220, command = self.gen_random)
        self.button_3.grid(row=4, column=0, padx=25, sticky="we")

        self.label_mode = customtkinter.CTkLabel(master=self.frame_left, text="Appearance Mode:")
        self.label_mode.grid(row=9, column=0, pady=0, padx=20, sticky="w")

        self.optionmenu_1 = customtkinter.CTkOptionMenu(master=self.frame_left,
                                                        values=["Light", "Dark", "System"],
                                                        command=self.change_appearance_mode)
        self.optionmenu_1.grid(row=10, column=0, pady=10, padx=20, columnspan=2, sticky="w")

        # ============ frame_right ============

        # configure grid layout (3x7)
        self.frame_right.rowconfigure((0, 1, 2, 3), weight=1)
        self.frame_right.rowconfigure(8, weight=10)
        self.frame_right.columnconfigure((0, 1), weight=1)
        self.frame_right.columnconfigure(2, weight=0)


        # ============ frame_info ============

        self.lab = customtkinter.CTkLabel(master = self.frame_right, image = self.img)
        self.lab.grid(row=0, column=0, columnspan = 3, rowspan=4)

        # ============ frame_right ============

        self.radio_var = tkinter.IntVar(value=0)

        self.label_radio_group = customtkinter.CTkLabel(master=self.frame_right,
                                                        text="Options")
        self.label_radio_group.grid(row=0, column=5, columnspan=1, pady=20, padx=10, sticky="")

        self.radio_button_1 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=0, text="Baby-step \ngiant-step", command = self.setpointcountmethod)
        self.radio_button_1.grid(row=1, column=5, pady=10, padx=20, sticky="n")

        self.radio_button_2 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=1, text="Schoof", command=self.setpointcountmethodSchoof)
        self.radio_button_2.grid(row=2, column=5, pady=10, padx=20, sticky="n")

        self.radio_button_3 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=2, text="SEA", command = self.setpointcountmethodSEA)
        self.radio_button_3.grid(row=3, column=5, pady=10, padx=20, sticky="n")

        self.radio_var1 = tkinter.IntVar(value=0)

        self.label_radio_group1 = customtkinter.CTkLabel(master=self.frame_right,
                                                        text="")
        self.label_radio_group1.grid(row=5, column=0, columnspan=1, pady=20, padx=10, sticky="")

        self.switch_1 = customtkinter.CTkSwitch(master=self.frame_right,
                                                text="Filter Atkin primes \n (recommended)",command = self.setfilter)
        self.switch_1.grid(row=4, column=5, columnspan=1, pady=10, padx=20, sticky="we")

        self.switch_2 = customtkinter.CTkSwitch(master=self.frame_right,
                                                text="Use Canonical modular\n polynomials (recommended)", command = self.setpoltype)
        self.switch_2.grid(row=5, column=5, columnspan=1, pady=10, padx=20, sticky="we")

        self.combobox_1 = customtkinter.CTkComboBox(master=self.frame_right,
                                                    values=["8", "16", "32", "64", "96", "128", "160", "192", "224", "256"], command = self.bitset)
        self.combobox_1.grid(row=4, column=1, columnspan=1, padx=20, sticky="we")


        self.sv = customtkinter.StringVar()
        self.sv.set("p")
        self.sv.trace("w", lambda name, index, mode=self.sv:self.getprimefromuser(self.sv))
        self.entry = customtkinter.CTkEntry(master=self.frame_right,
                                            width=120,
                                            placeholder_text=self.sv,textvariable=self.sv)

        self.sv2 = customtkinter.StringVar()
        self.sv2.trace("w", lambda name, index, mode=self.sv2: self.getAfromuser(self.sv2))
        self.entry.grid(row=5, column=0, columnspan=2,  padx=20, sticky="we")
        self.entry1 = customtkinter.CTkEntry(master=self.frame_right,
                                            width=120,
                                             textvariable = self.sv2,placeholder_text="A")
        self.entry1.grid(row=6, column=0, columnspan=2, pady = 5, padx=20, sticky="we")

        self.sv3 = customtkinter.StringVar()
        self.sv3.trace("w", lambda name, index, mode=self.sv3: self.getBfromuser(self.sv3))

        self.entry2 = customtkinter.CTkEntry(master=self.frame_right,
                                             width=120,
                                             placeholder_text="B", textvariable=self.sv3)
        self.entry2.grid(row=7, column=0, columnspan=2, pady = 5,padx=20, sticky="we")
        self.button_5 = customtkinter.CTkButton(master=self.frame_right,
                                                text="Calculate",
                                                border_width=2,  # <- custom border_width
                                                fg_color=None,  # <- no fg_color
                                                command=self.button_event)
        self.button_5.grid(row=7, column=5, columnspan=1, rowspan = 1, padx=20, sticky="we")

        self.output = customtkinter.CTkEntry(master=self.frame_right,placeholder_text="", width = 10000)


        self.entry1.insert(0, "A")
        self.entry2.insert(0, "B")
        self.output.grid(row=8, column=0, columnspan=6, rowspan = 2, padx=20,pady = 5, sticky="we")

        self.optionmenu_1.set("Dark")

        self.combobox_1.set("Bits")
        self.radio_button_3.select()
        self.clearbutton = customtkinter.CTkButton(master=self.frame_right,
                                                text="Clear",
                                                border_width=2,  # <- custom border_width
                                                fg_color=None,  # <- no fg_color
                                                command=self.clear)
        self.clearbutton.grid(row=6, column=5, columnspan=1, rowspan = 1, padx=20, sticky="we")
        self.switch_2.select()
        self.switch_1.select()
        self.output.configure(state="disabled")
        self.button_5.configure(state ="disabled")
        self.button_3.configure(state = "disabled")



    def gen_random(self):
        if self.combobox_1.get()== "Bits":
            self.output.configure(state="normal")
            self.output.delete(0, "end")
            self.output.configure(state="disabled")
            if not gmpy2.mpz(self.entry.textvariable.get()).is_prime():
                self.output.configure(state="normal")
                self.output.insert(0, "Error; p needs to be a prime")
                self.output.configure(state="disabled")
            else:
                self.userset = False
                self.prime = gmpy2.mpz(self.entry.textvariable.get())
                self.A, self.B = random.randint(1, self.prime - 1), random.randint(1, self.prime - 1)
                self.entry.textvariable.set(str(self.prime))
                self.entry1.textvariable.set(str(self.A))
                self.entry2.textvariable.set(str(self.B))

        else:
         prim = primefinder.randomnbitprime(int(self.primesize))
         self.prime = prim
         self.A, self.B = random.randint(1, self.prime-1),random.randint(1, self.prime-1)
         self.userset = False
         self.entry.textvariable.set(str(prim))

         self.entry1.textvariable.set(str(self.A))
         self.entry2.textvariable.set(str(self.B))

    def button_event(self):

        self.output.configure(state="normal")
        self.output.delete(0, "end")
        self.output.configure(state="disabled")
        if not gmpy2.mpz(self.entry.textvariable.get()).is_prime():

            self.output.configure(state = "normal")
            self.output.insert(0,"Error; p needs to be a prime")
            self.output.configure(state = "disabled")
        else:
            self.prime = gmpy2.mpz(self.entry.textvariable.get())
            self.A =  gmpy2.mpz(self.entry1.textvariable.get())
            self.B = gmpy2.mpz(self.entry2.textvariable.get())
            E = TestV2.EllipticCurve(self.prime, self.A, self.B)
            self.output.configure(state="normal")
            self.output.insert(0, "This may take a moment...")
            self.output.configure(state="disabled")


            self.button_5.configure(state = "disabled")
            self.button_3.configure(state = "disabled")
            self.entry.configure(state = "disabled")
            self.entry1.configure(state="disabled")
            self.entry2.configure(state="disabled")
            if self.method == 2:
              self.order = E.SEA(self.canonical, self.filter)
            if self.method == 1:
                self.order = E.Schoof()
            if self.method ==0:
                self.order == E.BabyStep()
            print(self.order)
            n = factorint(int(self.order[10:]))
            fac = ""
            for key in n:
                for i in range(n[key]):
                    fac = fac+str(key)+"*"
            fac = fac[:-1]
            divs =[]
            #for prime, multiplicity in n.items():
             #   divs += [prime] * multiplicity
           # Poi = E.findpoint(divs)
            self.output.configure(state="normal")
            self.output.delete(0, "end")
            self.output.configure(state="disabled")
            self.output.configure(state="normal")
            self.output.insert(0, self.order+"  =   "+fac)
            self.output.configure(state="disabled")
            self.entry.configure(state="normal")
            self.entry1.configure(state="normal")
            self.entry2.configure(state="normal")
    def bitset(self, bits):

        self.primesize = bits
        self.entry1.textvariable.set("A")
        self.entry2.textvariable.set("B")
        self.entry.textvariable.set("p")
        self.button_5.configure(state = "disabled")
        self.output.configure(state="normal")
        self.output.delete(0, "end")
        self.output.configure(state="disabled")
        self.button_3.configure(state="normal")
    def setfilter(self):

        self.filter = not self.filter
    def setpointcountmethodSEA(self):
        self.method = 2
        self.switch_1.configure(state = "normal")
        self.switch_2.configure(state="normal")
        self.clear()
    def setpointcountmethodSchoof(self):
        self.method = 1
        self.switch_1.configure(state="disabled")
        self.switch_2.configure(state="disabled")
        self.clear()
    def setpointcountmethod(self):
        self.method = 0
        self.switch_1.configure(state="disabled")
        self.switch_2.configure(state="disabled")
        self.clear()
    def button_event1(self):
        print("nothing here yet")
    def setpoltype(self):

        self.canonical = not self.canonical

    def change_appearance_mode(self, new_appearance_mode):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def on_closing(self, event=0):
        self.destroy()
    def getprimefromuser(self, other):
        #if other.get()!="p":
        #  self.combobox_1.set("Bits")
        if primefinder.checktype(other.get()) and self.userset == True:
            self.combobox_1.set("Bits")
        if other is None or other.get() == "" or not primefinder.checktype(other.get()):

            self.button_5.configure(state="disabled")
            self.button_3.configure(state="disabled")
        else:

         self.button_3.configure(state="normal")
         self.prime = other.get()

         if self.entry1.textvariable is not None and self.entry2.textvariable is not None:
            if self.entry2.textvariable.get() not in ["", "B"] and self.entry1.textvariable.get() not in ["", "A"]:
                self.button_5.configure(state="normal")
                self.userset = True

    def getAfromuser(self, other):
        #if other.get()!= "A":
            #self.combobox_1.set("Bits")
        self.A = other.get()

        if other is None or other.get() == "" or not primefinder.checktype(other.get()):

            self.button_5.configure(state = "disabled")
        else:
            if self.entry2.textvariable is not None and self.entry.textvariable is not None:
                 if self.entry2.textvariable.get() not in [ "", "B"] and self.entry.textvariable.get()not in ["", "p"]:
                     self.button_5.configure(state = "normal")

    def getBfromuser(self, other):
        #if other.get()!= "B":
       #  self.combobox_1.set("Bits")
        self.B = other.get()
        if other is None or other.get() == "" or not primefinder.checktype(other.get()):

            self.button_5.configure(state="disabled")
        else:
            if self.entry1.textvariable is not None and self.entry.textvariable is not None:
                 if self.entry1.textvariable.get() not in ["", "A"] and self.entry.textvariable.get() not in ["", "p"]:
                     self.button_5.configure(state = "normal")

    def clear(self):
        self.A,self.B, self.prime, self.primesize = 0,0,0,0
        self.button_5.configure(state = "disabled")
        self.button_3.configure(state = "disabled")
        self.combobox_1.set("Bits")
        self.entry.textvariable.set("p")
        self.entry1.textvariable.set("A")
        self.entry2.textvariable.set("B")
        self.output.configure(state="normal")
        self.output.delete(0, "end")
        self.output.configure(state="disabled")
        self.userset = True
        self.order = 0
        #self.method = 2
        #self.radio_button_3.select()


if __name__ == "__main__":
    app = App()
    app.mainloop()
