
simply models t->bW->qq

To run:
root -q run.C

In the output root tree, you will get
q1p4, q2p4, Wp4, bp4, topp4, dRqq, dRq1b, dRq2b branches

Examples:
$ root data.root
root [0] Events->Draw("dRqq")
root [0] Events->Draw("topp4.E()")
