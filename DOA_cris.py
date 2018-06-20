import numpy
from gnuradio import gr

class doa_cris(gr.sync_block):
    """
    docstring for block oreja
    """
    def __init__(self):
        gr.sync_block.__init__(self,
            name="oreja",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])
        self.muestras = []
        self.M = 2
        self.d = 0.26
        self.f = 915e6
        self.Rxx = []


    def work(self, input_items, output_items):
        output_items[0][:] = 0*input_items[0]
        xx = numpy.matrix((input_items[0],input_items[0])) 
        self.Rxx = numpy.matmul(numpy.asmatrix(xx),numpy.asmatrix(xx).H)
        #print "el xx es " + str(numpy.asmatrix(xx))
        #print "el xx es " + str(numpy.asmatrix(xx).H)
        #print  "la matriz " + str(self.Rxx)
        #print "el tamano de la entrada es" + str(len(input_items[0]))
        #print "el tamano de Rxx es " + str(len(self.Rxx))
        Rxx = self.Rxx
        theta= range(181)
        th=numpy.ones(180)
        P_Bartlett=numpy.ones(180)
        P_Capon=numpy.ones(180)
        P_MUSIC= numpy.ones(180)
        for k in range(180):
            th[k]=theta[k]*numpy.pi/180
            a=1
            for jj in range(2,self.M+1):
                a =numpy.matrix([a,numpy.exp(1j*jj*numpy.pi*numpy.sin(th[k]))])
            P_Bartlett[k]=numpy.real(numpy.matmul(numpy.matmul(a,Rxx),a.T))
        output_items[0][0:180] = P_Bartlett
        self.muestras = input_items[0]
        return len(output_items[0])