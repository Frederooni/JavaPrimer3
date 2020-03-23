/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Converted to Java by Fred Long */

package fred.bio.primer3;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import static java.lang.Math.*;
import static java.util.Arrays.binarySearch;

public class Thal {

    /** Set this to true to output extra debugging information. */
    public boolean DEBUG = false;

    /**
     * The maximum length of _one_ of the two sequences being aligned in a thermodynamic alignment. In other words, the
     * length of one sequence must be <= THAL_MAX_ALIGN, but the other sequence can be longer. The rationale behind this
     * value (60) is that this is the maxium reasonable length for nearest neighbor models. It is the maxium length at
     * which we can restrict our model to only two states of melting: fully intact duplex or completely dissociated
     * single strands.
     */
    private static final int THAL_MAX_ALIGN = 60;

    /**
     * The maxium length of the other sequence in a thermodynamic alignment. This value can be increased, though
     * alignments against very long sequences will be quite slow. As of 2012-05-18, we only potentially see sequences
     * longer this when checking for mispriming in the template ('max_template_mispriming') in libprimer3.c, which is
     * really designed to find sites of ectopic primer very close (a few kilobases) from the location of the cadidate
     * primer.
     */
    private static final int THAL_MAX_SEQ = 10000;

    public static final int THAL_ANY = 1;
    public static final int THAL_END1 = 2;
    public static final int THAL_END2 = 3;
    public static final int THAL_HAIRPIN = 4;

    /**
     * Structure for passing arguments to THermodynamic ALignment calculation
     */
    public static class ThalArgs {
        /** Output mode */
        public int mode = OUTPUT_MODE_FAST;
        /**
         * Alignment type:
         * THAL_ANY (default), THAL_END1, THAL_END2, or THAL_HAIRPIN
         */
        public int type = THAL_ANY;
        /**
         * maximum size of loop to consider; longer than 30 bp are not allowed
         */
        public int maxLoop = MAX_LOOP;
        /**
         * mM concentration of monovalent cations (default 50 mM)
         */
        public double mv = 50;
        /**
         * mM concentration of divalent cations (default 1.5 mM)
         */
        public double dv = 1.5;
        /**
         * mM concentration of dNTP-s (default 0.6 mM)
         */
        public double dntp = 0.6;
        /**
         * nM concentration of oligonucleotides (default 50 nM)
         */
        public double dna_conc = 50;
        /**
         * Kelvin temperature from which hairpin structures will be calculated (default 310.15 K, 37 C)
         */
        public double temp = TEMP_KELVIN;
        /** calculate dimer structure (no longer used) */
        public boolean dimer = true;
    }

    /* Structure for receiving results from the thermodynamic alignment calculation */
    public static class ThalResults {
        /** Informational message */
        String msg;
        /** Melting temperature in Celsius */
        public double temp;
        /** Delta G (\u0394G) */
        double delta_G;
        /** Delta H (\u0394H) */
        double delta_H;
        /** Delta S (\u0394S) */
        double delta_S;
        int align_end_1;
        int align_end_2;
        /** Secondary structure */
        String sec_struct;
    }

    /**
     * Don't compute a secondary structure string in {@link ThalResults#sec_struct}.
     */
    public static final int OUTPUT_MODE_FAST = 0;
    /**
     * Return a text-based secondary structure string in {@link ThalResults#sec_struct}.
     */
    public static final int OUTPUT_MODE_TEXT = 1;
    /**
     * Return a graphical secondary structure string in {@link ThalResults#sec_struct}.
     */
    public static final int OUTPUT_MODE_GRAPHICAL = 4;

    /**  minimum size of hairpin loop */
    private static final int MIN_HRPN_LOOP = 3;

    /* BEGIN CONSTANTS */

    private static final double SMALL_NON_ZERO = 0.000001; // 1e-6
    private static final double _INFINITY = 1.0 / 0.0;
    private static final double THAL_ERROR_SCORE = -_INFINITY;
    /** cal/Kmol */
    private static final double R = 1.9872;
    /** Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
    private static final double ILAS = (-300 / 310.15);
    /** Internal Loop EntHalpy Asymmetry correction */
    private static final double ILAH = 0.0;
    /** AT penalty */
    private static final double AT_H = 2200.0;
    /** AT penalty */
    private static final double AT_S = 6.9;
    /** to filter out non-existing entropies */
    private static final double MinEntropyCutoff = -2500.0;
    /** initiation */
    private static final double MinEntropy = -3224.0;
    /** structures w higher G are considered to be unstabile */
    private static final double G2 = 0.0;
    private static final double ABSOLUTE_ZERO = 273.15;
    /** 310.15 Kelvin is 37 Celsius */
    private static final double TEMP_KELVIN = 310.15;
    /** the maximum size of loop that can be calculated; for larger loops formula must be implemented */
    private static final int MAX_LOOP = 30;
    private static final char BASES[] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
    /** default thermodynamic parameters */
    private static ThalParameters defaultParameters = new ThalParameters();

    /** matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
    private static final int BPI[][] = {
            {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
            {0, 0, 1, 0, 0},
            {0, 1, 0, 0, 0},
            {1, 0, 0, 0, 0},
            {0, 0, 0, 0, 0}};

    /* END OF CONSTANTS */

    /* BEGIN STRUCTs */

    private static class loop {
        byte loop[];
        double value;
        loop(int size) {
            this.loop = new byte[size];
        }
        loop(byte[] loop) {
            this.loop = loop;
        }
    }

    private static class triloop extends loop {
        triloop() {
            super(5);
        }
    }

    private static class tetraloop extends loop {
        tetraloop() {
            super(6);
        }
    }

    /** structure for tracebacku - unimolecular str */
    private static class tracer {
        int i;
        int j;
        int mtrx; /* [0 1] entropyDPT/enthalpyDPT*/
        tracer next;
    }

    /* END STRUCTs */

    private static ThreadLocal<Thal> thalHolder = new ThreadLocal<>();

    /** AT penalty */
    private double atPenaltyS[][] = new double[5][5];
    /** AT penalty */
    private double atPenaltyH[][] = new double[5][5];
    /* entropies of most stable hairpin terminal bp at 5' end */
    private double send5[];
    /* enthalpies of most stable hairpin terminal bp at 5' end */
    private double hend5[];
    /* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
    /** initiation enthalpy; for duplex 200, for unimolecular structure 0 */
    private double dplx_init_H;
    /** initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
    private double dplx_init_S;
    /** value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
    private double saltCorrection;
    /** universal gas constant multiplied w DNA conc - for melting temperature */
    private double RC;
    /** var that helps to find str w highest melting temperature */
    private double SHleft;
    /** starting position of most stable str */
    private int bestI, bestJ;
    /* dynamic programming matrix for values of enthalpy */
    private double enthalpyDPT[][];
    /* dynamic programming matrix for values of entropy */
    private double entropyDPT[][];
    /** Input oligo sequences */
    private String oligo1, oligo2;
    /** same as oligo1 and oligo2 but converted to numbers */
    private byte[] numSeq1, numSeq2;
    /** length of sequense 1 and 2 */
    private int len1, len2, len3;
    /** thermodynamic paramteres for 3' dangling ends */
    private double dangleEntropies3[][][] = new double[5][5][5];
    /** ther params for 3' dangling ends */
    private double dangleEnthalpies3[][][] = new double[5][5][5];
    /** ther params for 5' dangling ends */
    private double dangleEntropies5[][][] = new double[5][5][5];
    /** ther params for 5' dangling ends */
    private double dangleEnthalpies5[][][] = new double[5][5][5];
    /** ther params for perfect match pairs */
    private double stackEntropies[][][][] = new double[5][5][5][5];
    /** ther params for perfect match pairs */
    private double stackEnthalpies[][][][] = new double[5][5][5][5];
    /** ther params for perfect match and internal mm */
    private double stackint2Entropies[][][][] = new double[5][5][5][5];
    /** ther params for perfect match and internal mm*/
    private double stackint2Enthalpies[][][][] = new double[5][5][5][5];
    /** interior loop params according to length of the loop */
    private double interiorLoopEntropies[] = new double[30];
    /** bulge loop params according to length of the loop */
    private double bulgeLoopEntropies[] = new double[30];
    /** hairpin loop params according to length of the loop */
    private double hairpinLoopEntropies[] = new double[30];
    /** same as interiorLoopEntropies but values of entropy */
    private double interiorLoopEnthalpies[] = new double[30];
    /** same as bulgeLoopEntropies but values of entropy */
    private double bulgeLoopEnthalpies[] = new double[30];
    /** same as hairpinLoopEntropies but values of entropy */
    private double hairpinLoopEnthalpies[] = new double[30];
    /** ther params for terminal mismatches */
    private double tstackEntropies[][][][] = new double[5][5][5][5];
    /** ther params for terminal mismatches */
    private double tstackEnthalpies[][][][] = new double[5][5][5][5];
    /** ther params for internal terminal mismatches */
    private double tstack2Entropies[][][][] = new double[5][5][5][5];
    /** ther params for internal terminal mismatches */
    private double tstack2Enthalpies[][][][] = new double[5][5][5][5];
    /** number of hairpin triloop penalties */
    private int numTriloops;
    /** ther penalties for given triloop seq-s */
    private triloop triloopEntropies[] = null;
    /** ther penalties for given triloop seq-s */
    private triloop triloopEnthalpies[] = null;
    /** number of hairpin tetraloop penalties */
    private int numTetraloops;
    /** ther penalties for given tetraloop seq-s */
    private tetraloop tetraloopEntropies[] = null;
    /** ther penalties for given tetraloop seq-s */
    private tetraloop tetraloopEnthalpies[] = null;
    /** the current set to thermodynamic parameters */
    private ThalParameters thal_parameters;

    /** Create Thal instance with the given thermodynamic parameters. */
    public Thal(ThalParameters tp) throws Exception {
        get_thermodynamic_values(tp);
    }

    public void get_thermodynamic_values(ThalParameters tp) throws Exception {
        this.thal_parameters = tp;
        /* Read the thermodynamic values (parameters) from the parameter files
           in the directory specified by 'path'.  Return 0 on success and -1
           on error. The thermodynamic values are stored in multiple static
           variables. */
        getStack(stackEntropies, stackEnthalpies, tp);
        verifyStackTable(stackEntropies, "entropy"); /* this is for code debugging */
        verifyStackTable(stackEnthalpies, "enthalpy"); /* this is for code debugging */
        getStackint2(stackint2Entropies, stackint2Enthalpies, tp);
        getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5, tp);
        getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
                interiorLoopEnthalpies, bulgeLoopEnthalpies, tp);
        getTstack(tstackEntropies, tstackEnthalpies, tp);
        getTstack2(tstack2Entropies, tstack2Enthalpies, tp);
        getTriloop(tp);
        getTetraloop(tp);
        /* getting the AT-penalties */
        tableStartATS(AT_S, atPenaltyS);
        tableStartATH(AT_H, atPenaltyH);
    }

    /** Create Thal instance with default thermodynamic parameters. */
    public Thal() throws Exception {
        this(defaultParameters);
    }

    /** Creates a thread-local instance of Thal for the current thread. */
    public static Thal getThreadLocalThal() throws Exception {
        Thal thal = thalHolder.get();
        if (thal == null) {
            thal = new Thal();
            thalHolder.set(thal);
        }
        return thal;
    }

    /** Creates a thread-local instance of Thal for the current thread. */
    public static Thal getThreadLocalThal(ThalParameters p) throws Exception {
        Thal thal = thalHolder.get();
        if (thal == null) {
            thal = new Thal(p);
            thalHolder.set(thal);
        }
        if (thal.thal_parameters != p) thal.get_thermodynamic_values(p);
        return thal;
    }

    /**
     * Computes the melting temperature, secondary structure, and other thermodynamic values for the two oligos.
     * If calculating hairpins, oligo_r should be the same as oligo_l.
     * <p>
     * This method is not thread-safe.  To call this method from multiple threads, create a separate
     * Thal instance for each thread, or use {@link #getThreadLocalThal()}, which will return a thread-local
     * instance of Thal for the current thread.  Usage example:
     * <pre>
     *     Thal.getThreadLocalThal.thal(oligo_l, oligo_r, args);
     * </pre>
     */
    public ThalResults thal(String oligo_f, String oligo_r, ThalArgs a) throws RuntimeException {
        CHECK_ERROR(null == a, "null args");
        ThalResults o = new ThalResults();
        int mode = a.mode;
        double[] SH;
        int i, j;
        int len_f, len_r;
        int k;
        int[] bp;
        String oligo2_rev = null;
        double mh, ms;
        double G1, bestG;

        send5 = hend5 = null;
        enthalpyDPT = entropyDPT = null;
        numSeq1 = numSeq2 = null;
        oligo1 = oligo2 = null;
        o.msg = "";
        o.temp = THAL_ERROR_SCORE;

        CHECK_ERROR(null == oligo_f, "null first sequence");
        CHECK_ERROR(null == oligo_r, "null second sequence");
        len_f = length_unsig_char(oligo_f);
        len_r = length_unsig_char(oligo_r);
        CHECK_ERROR(len_f == 0, "empty first sequence");
        CHECK_ERROR(len_r == 0, "empty second sequence");

       /*CHECK_ERROR(1==len_f, "Length 1 first sequence");
       CHECK_ERROR(1==len_r, "Length 1 second sequence"); */
       /* The following error messages will be seen by end users and will
          not be easy to understand. */
        CHECK_ERROR((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN),
                "Both sequences longer than " + THAL_MAX_ALIGN + " for thermodynamic alignment");
        CHECK_ERROR((len_f > THAL_MAX_SEQ), LONG_SEQ_ERR_STR(THAL_MAX_SEQ) + " (1)");
        CHECK_ERROR((len_r > THAL_MAX_SEQ), LONG_SEQ_ERR_STR(THAL_MAX_SEQ) + " (2)");

        CHECK_ERROR(a.type < THAL_ANY || a.type > THAL_HAIRPIN, "Illegal alignment type");
        o.align_end_1 = -1;
        o.align_end_2 = -1;
        if (a.type != THAL_END2) {
            oligo1 = oligo_f;
            oligo2 = oligo_r;
        } else {
            oligo1 = oligo_r;
            oligo2 = oligo_f;
        }
        /* INIT values for unimolecular and bimolecular structures */
        if (a.type == THAL_HAIRPIN) { /* unimolecular folding */
            len2 = oligo2.length();
            len3 = len2 - 1;
            dplx_init_H = 0.0;
            dplx_init_S = -0.00000000001;
            RC = 0;
        } else if (a.type != THAL_HAIRPIN) {
            /* hybridization of two oligos */
            dplx_init_H = 200;
            dplx_init_S = -5.7;
            if (symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
                RC = R * log(a.dna_conc / 1000000000.0);
            } else {
                RC = R * log(a.dna_conc / 4000000000.0);
            }
            /* REVERSE oligo2, so it goes to dpt 3'.5' direction */
            if (a.type != THAL_END2) {
                oligo2_rev = reverse(oligo_r);
            } else {
                oligo2_rev = reverse(oligo_f);
            }
            oligo2 = oligo2_rev;
        }
        len1 = oligo1.length();
        len2 = oligo2.length();
        /* convert nucleotides to numbers */
        numSeq1 = new byte[len1 + 2];
        numSeq2 = new byte[len2 + 2];

        /* Calc part of the salt correction */
        /* salt correction for entropy, must be multiplied with N, which is the total number of phosphates
           in the duplex divided by 2; 8bp dplx N=7 */
        saltCorrection = saltCorrectS(a.mv, a.dv, a.dntp);

        enthalpyDPT = new double[len1 + 1][len2 + 1]; /* dyn. programming table for dS and dH */
        entropyDPT = new double[len1 + 1][len2 + 1]; /* enthalpyDPT is 3D array represented as 1D array */
        if (a.type == THAL_HAIRPIN) { /* monomer */
            /* terminal basepairs */
            send5 = new double[len1 + 1];
            hend5 = new double[len1 + 1];
        }
        oligo1 = oligo1.toUpperCase();
        oligo2 = oligo2.toUpperCase();
        for (i = 1; i <= len1; ++i) numSeq1[i] = str2int(oligo1.charAt(i - 1));
        for (i = 1; i <= len2; ++i) numSeq2[i] = str2int(oligo2.charAt(i - 1));
        numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4; /* mark as N-s */
        if (a.type == THAL_HAIRPIN) { /* calculate structure of monomer */
            initMatrix2();
            fillMatrix2(a.maxLoop, o);
            calc_terminal_bp(a.temp);
            mh = hend5[len1];
            ms = send5[len1];
            o.align_end_1 = (int) mh;
            o.align_end_2 = (int) ms;
            bp = new int[len1];
            for (k = 0; k < len1; ++k) bp[k] = 0;
            if (isFinite(mh)) {
                tracebacku(bp, a.maxLoop, o);
                /* traceback for unimolecular structure */
                o.sec_struct=drawHairpin(bp, mh, ms, mode,a.temp, o); /* if mode=THL_FAST or THL_DEBUG_F then return after printing basic therm data */
            } else if ((mode != OUTPUT_MODE_FAST) && (mode != OUTPUT_MODE_GRAPHICAL)) {
                System.err.print("No secondary structure could be calculated\n"); // TODO: Throw an exception?
            }

            if (o.temp == -_INFINITY && (empty_string(o.msg))) o.temp = 0.0;
            bp = null;
            enthalpyDPT = null;
            entropyDPT = null;
            numSeq1 = null;
            numSeq2 = null;
            send5 = null;
            hend5 = null;
            oligo1 = null;
            oligo2 = null;
            return o;
        } else if (a.type != THAL_HAIRPIN) { /* Hybridization of two moleculs */
            len3 = len2;
            initMatrix();
            fillMatrix(a.maxLoop, o);
            SH = new double[2];
            /* calculate terminal basepairs */
            bestI = bestJ = 0;
            G1 = bestG = _INFINITY;
            if (a.type == THAL_ANY)
                for (i = 1; i <= len1; i++) {
                    for (j = 1; j <= len2; j++) {
                        RSH(i, j, SH);
                        SH[0] = SH[0] + SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
                        SH[1] = SH[1] + SMALL_NON_ZERO;
                        G1 = (enthalpyDPT[i][j] + SH[1] + dplx_init_H) - TEMP_KELVIN * (entropyDPT[i][j] + SH[0] + dplx_init_S);
                        if (G1 < bestG) {
                            bestG = G1;
                            bestI = i;
                            bestJ = j;
                        }
                    }
                }
            int[] ps1 = new int[len1];
            int[] ps2 = new int[len2];
            for (i = 0; i < len1; ++i) ps1[i] = 0;
            for (j = 0; j < len2; ++j) ps2[j] = 0;
            if (a.type == THAL_END1 || a.type == THAL_END2) {
                /* THAL_END1 */
                bestI = bestJ = 0;
                bestI = len1;
                i = len1;
                G1 = bestG = _INFINITY;
                for (j = 1; j <= len2; ++j) {
                    RSH(i, j, SH);
                     /* this adding is done for compiler, optimization -O2 vs -O0,
                        that compiler could understand that SH is changed in this cycle */
                    SH[0] = SH[0] + SMALL_NON_ZERO;
                    SH[1] = SH[1] + SMALL_NON_ZERO;
                    G1 = (enthalpyDPT[i][j] + SH[1] + dplx_init_H) - TEMP_KELVIN * (entropyDPT[i][j] + SH[0] + dplx_init_S);
                    if (G1 < bestG) {
                        bestG = G1;
                        bestJ = j;
                    }
                }
            }
            if (!isFinite(bestG)) bestI = bestJ = 1;
            double dH, dS;
            RSH(bestI, bestJ, SH);
            dH = enthalpyDPT[bestI][bestJ] + SH[1] + dplx_init_H;
            dS = (entropyDPT[bestI][bestJ] + SH[0] + dplx_init_S);
            /* tracebacking */
            for (i = 0; i < len1; ++i) ps1[i] = 0;
            for (j = 0; j < len2; ++j) ps2[j] = 0;
            if (isFinite(enthalpyDPT[bestI][bestJ])) {
                traceback(bestI, bestJ, RC, ps1, ps2, a.maxLoop, o);
                o.sec_struct = drawDimer(ps1, ps2, SHleft, dH, dS, mode, a.temp, o);
                o.align_end_1 = bestI;
                o.align_end_2 = bestJ;
            } else {
                o.temp = 0.0;
                /* fputs("No secondary structure could be calculated\n",stderr); */
            }
            ps1 = null;
            ps2 = null;
            SH = null;
            oligo2_rev = null;
            enthalpyDPT = null;
            entropyDPT = null;
            numSeq1 = null;
            numSeq2 = null;
            oligo1 = null;
            return o;
        }
        return o;
    }
    /*** END thal() ***/

    /**
     * converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever
     */
    private byte
    str2int(char c) {
        switch (c) {
            case 'A':
            case '0':
                return 0;
            case 'C':
            case '1':
                return 1;
            case 'G':
            case '2':
                return 2;
            case 'T':
            case '3':
                return 3;
        }
        return 4;
    }


    private int max5(double a, double b, double c, double d, double e) {
        if (a > b && a > c && a > d && a > e) return 1;
        else if (b > c && b > d && b > e) return 2;
        else if (c > d && c > e) return 3;
        else if (d > e) return 4;
        else return 5;
    }

    /**
     * to add elements to struct
     */
    private tracer push(tracer stack, int i, int j, int mtrx) {
        tracer new_top = new tracer();
        new_top.i = i;
        new_top.j = j;
        new_top.mtrx = mtrx;
        new_top.next = stack;
        return new_top;
    }

    private String reverse(String s) {
        StringBuilder sb = new StringBuilder(s);
        return sb.reverse().toString();
    }

    /**
     * part of calculating salt correction for Tm by SantaLucia et al
     */
    private double
    saltCorrectS(double mv, double dv, double dntp) {
        if (dv <= 0) dntp = dv;
        return 0.368 * ((log((mv + 120 * (sqrt(max(0.0, dv - dntp)))) / 1000)));
    }

    /**
     * These functions are needed as "inf" cannot be read on Windows directly
     */
    private double
    readDouble(BufferedReader str) throws Exception {
        return string_to_double(str.readLine().trim());
    }

    private double string_to_double(String str) {
        if (str.startsWith("inf")) return _INFINITY;
        else return Double.parseDouble(str);
    }

    /**
     * Reads a line containing 4 doubles, which can be specified as "inf".
     */
    private void
    readLoop(BufferedReader str, int k, double[] hairpinLoop, double[] interiorLoop, double[] bulgeLoop) throws Exception {
        String line = str.readLine().trim();
        String[] numbers = line.split("\\s+");
        int index = Integer.parseInt(numbers[0]);
        assert index == k + 1;
        interiorLoop[k] = string_to_double(numbers[1]);
        bulgeLoop[k] = string_to_double(numbers[2]);
        hairpinLoop[k] = string_to_double(numbers[3]);
    }

    /**
     * Reads a line containing a short string and a double, used for reading a triloop or tetraloop.
     */
    private boolean
    readTLoop(BufferedReader str, StringBuilder sb, double[] v) throws IOException {
        String line = str.readLine();
        if (line == null) return false;
        String[] words = line.trim().split("\\s+");
        sb.delete(0, sb.length());
        sb.append(words[0]);
        v[0] = string_to_double(words[1]);
        return true;
    }

    private void
    getStack(double stackEntropies[][][][], double stackEnthalpies[][][][], final ThalParameters tp) throws Exception {
        int i, j, ii, jj;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.stack_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.stack_dh));
        for (i = 0; i < 5; ++i) {
            for (ii = 0; ii < 5; ++ii) {
                for (j = 0; j < 5; ++j) {
                    for (jj = 0; jj < 5; ++jj) {
                        if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                            stackEntropies[i][ii][j][jj] = -1.0;
                            stackEnthalpies[i][ii][j][jj] = _INFINITY;
                        } else {
                            stackEntropies[i][ii][j][jj] = readDouble(pt_ds);
                            stackEnthalpies[i][ii][j][jj] = readDouble(pt_dh);
                            if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
                                stackEntropies[i][ii][j][jj] = -1.0;
                                stackEnthalpies[i][ii][j][jj] = _INFINITY;
                            }
                        }
                    }
                }
            }
        }
    }

    private void
    getStackint2(double stackint2Entropies[][][][], double stackint2Enthalpies[][][][], final ThalParameters tp) throws Exception {
        int i, j, ii, jj;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.stackmm_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.stackmm_dh));
        for (i = 0; i < 5; ++i) {
            for (ii = 0; ii < 5; ++ii) {
                for (j = 0; j < 5; ++j) {
                    for (jj = 0; jj < 5; ++jj) {
                        if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                            stackint2Entropies[i][ii][j][jj] = -1.0;
                            stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                        } else {
                            stackint2Entropies[i][ii][j][jj] = readDouble(pt_ds);
                            stackint2Enthalpies[i][ii][j][jj] = readDouble(pt_dh);
                            if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                                stackint2Entropies[i][ii][j][jj] = -1.0;
                                stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                            }
                        }
                    }
                }
            }
        }
    }

    private void verifyStackTable(double stack[][][][], String type) {
        boolean error = false;
        int i, j, ii, jj;
        for (i = 0; i < 4; ++i)
            for (j = 0; j < 4; ++j)
                for (ii = 0; ii < 4; ++ii)
                    for (jj = 0; jj < 4; ++jj)
                        if (stack[i][j][ii][jj] != stack[jj][ii][j][i]) {
                            System.err.printf("Warning: symmetrical stacks _are_ _not_ equal: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n",
                                    BASES[i], BASES[j], BASES[ii], BASES[jj], type, stack[i][j][ii][jj], BASES[jj],
                                    BASES[ii], BASES[j], BASES[i], type, stack[jj][ii][j][i]);
                            error = true;
                        }
        CHECK_ERROR(error, "Non-symmetry found in %s matrix " + type);
    }

    private void
    getDangle(double dangleEntropies3[][][], double dangleEnthalpies3[][][], double dangleEntropies5[][][],
              double dangleEnthalpies5[][][], final ThalParameters tp) throws Exception {
        int i, j, k;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.dangle_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.dangle_dh));
        for (i = 0; i < 5; ++i)
            for (j = 0; j < 5; ++j)
                for (k = 0; k < 5; ++k) {
                    if (i == 4 || j == 4) {
                        dangleEntropies3[i][k][j] = -1.0;
                        dangleEnthalpies3[i][k][j] = _INFINITY;
                    } else if (k == 4) {
                        dangleEntropies3[i][k][j] = -1.0;
                        dangleEnthalpies3[i][k][j] = _INFINITY;
                    } else {
                        dangleEntropies3[i][k][j] = readDouble(pt_ds);
                        dangleEnthalpies3[i][k][j] = readDouble(pt_dh);
                        if (!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
                            dangleEntropies3[i][k][j] = -1.0;
                            dangleEnthalpies3[i][k][j] = _INFINITY;
                        }
                    }
                }

        for (i = 0; i < 5; ++i)
            for (j = 0; j < 5; ++j)
                for (k = 0; k < 5; ++k) {
                    if (i == 4 || j == 4) {
                        dangleEntropies5[i][j][k] = -1.0;
                        dangleEnthalpies5[i][j][k] = _INFINITY;
                    } else if (k == 4) {
                        dangleEntropies5[i][j][k] = -1.0;
                        dangleEnthalpies5[i][j][k] = _INFINITY;
                    } else {
                        dangleEntropies5[i][j][k] = readDouble(pt_ds);
                        dangleEnthalpies5[i][j][k] = readDouble(pt_dh);
                        if (!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
                            dangleEntropies5[i][j][k] = -1.0;
                            dangleEnthalpies5[i][j][k] = _INFINITY;
                        }
                    }
                }
    }

    private void
    getLoop(double hairpinLoopEntropies[], double interiorLoopEntropies[], double bulgeLoopEntropies[],
            double hairpinLoopEnthalpies[], double interiorLoopEnthalpies[], double bulgeLoopEnthalpies[],
            final ThalParameters tp) throws Exception {
        int k;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.loops_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.loops_dh));
        for (k = 0; k < 30; ++k) {
            readLoop(pt_ds, k, hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies);
            readLoop(pt_dh, k, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies);
        }
    }

    private void
    getTstack(double tstackEntropies[][][][], double tstackEnthalpies[][][][], final ThalParameters tp) throws Exception {
        int i1, j1, i2, j2;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.tstack_tm_inf_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.tstack_dh));
        for (i1 = 0; i1 < 5; ++i1)
            for (i2 = 0; i2 < 5; ++i2)
                for (j1 = 0; j1 < 5; ++j1)
                    for (j2 = 0; j2 < 5; ++j2)
                        if (i1 == 4 || j1 == 4) {
                            tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
                            tstackEntropies[i1][i2][j1][j2] = -1.0;
                        } else if (i2 == 4 || j2 == 4) {
                            tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
                            tstackEnthalpies[i1][i2][j1][j2] = 0.0;
                        } else {
                            tstackEntropies[i1][i2][j1][j2] = readDouble(pt_ds);
                            tstackEnthalpies[i1][i2][j1][j2] = readDouble(pt_dh);
                            if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
                                tstackEntropies[i1][i2][j1][j2] = -1.0;
                                tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
                            }
                        }
    }

    private void
    getTstack2(double tstack2Entropies[][][][], double tstack2Enthalpies[][][][], final ThalParameters tp) throws Exception {

        int i1, j1, i2, j2;
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.tstack2_ds));
        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.tstack2_dh));
        for (i1 = 0; i1 < 5; ++i1)
            for (i2 = 0; i2 < 5; ++i2)
                for (j1 = 0; j1 < 5; ++j1)
                    for (j2 = 0; j2 < 5; ++j2)
                        if (i1 == 4 || j1 == 4) {
                            tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
                            tstack2Entropies[i1][i2][j1][j2] = -1.0;
                        } else if (i2 == 4 || j2 == 4) {
                            tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
                            tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
                        } else {
                            tstack2Entropies[i1][i2][j1][j2] = readDouble(pt_ds);
                            tstack2Enthalpies[i1][i2][j1][j2] = readDouble(pt_dh);
                            if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
                                tstack2Entropies[i1][i2][j1][j2] = -1.0;
                                tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
                            }
                        }
    }

    private void
    getTriloop(final ThalParameters tp) throws IOException {
        double[] value = new double[1];
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.triloop_ds));
        numTriloops = 0;
        ArrayList<triloop> triloops = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        while (readTLoop(pt_ds, sb, value)) {
            triloop loop = new triloop();
            for (int i = 0; i < 5; ++i) loop.loop[i] = str2int(sb.charAt(i));
            loop.value = value[0];
            triloops.add(loop);
            ++numTriloops;
        }
        triloopEntropies = triloops.toArray(new triloop[0]);

        BufferedReader pt_dh = new BufferedReader(new StringReader(tp.triloop_dh));
        numTriloops = 0;
        triloops.clear();
        while (readTLoop(pt_dh, sb, value)) {
            triloop loop = new triloop();
            for (int i = 0; i < 5; ++i) loop.loop[i] = str2int(sb.charAt(i));
            loop.value = value[0];
            triloops.add(loop);
            ++numTriloops;
        }
        triloopEnthalpies = triloops.toArray(new triloop[0]);
    }

    private void
    getTetraloop(final ThalParameters tp) throws IOException {
        double[] value = new double[1];
        BufferedReader pt_ds = new BufferedReader(new StringReader(tp.tetraloop_ds));
        numTetraloops = 0;
        ArrayList<tetraloop> tetraloops = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        while (readTLoop(pt_ds, sb, value)) {
            tetraloop loop = new tetraloop();
            for (int i = 0; i < 6; ++i) loop.loop[i] = str2int(sb.charAt(i));
            loop.value = value[0];
            tetraloops.add(loop);
            ++numTetraloops;
        }
        tetraloopEntropies = tetraloops.toArray(new tetraloop[0]);

        numTetraloops = 0;
        tetraloops.clear();
        while (readTLoop(pt_ds, sb, value)) {
            tetraloop loop = new tetraloop();
            for (int i = 0; i < 6; ++i) loop.loop[i] = str2int(sb.charAt(i));
            loop.value = value[0];
            tetraloops.add(loop);
            ++numTetraloops;
        }
        tetraloopEnthalpies = tetraloops.toArray(new tetraloop[0]);
    }

    /**
     * creates table of entropy values for nucleotides to which AT-penlty must be applied
     */
    private void
    tableStartATS(double atp_value, double atPenaltyS[][]) {
        int i, j;
        for (i = 0; i < 5; ++i)
            for (j = 0; j < 5; ++j)
                atPenaltyS[i][j] = 0.00000000001;
        atPenaltyS[0][3] = atPenaltyS[3][0] = atp_value;
    }


    private void
    tableStartATH(double atp_value, double atPenaltyH[][]) {
        int i, j;
        for (i = 0; i < 5; ++i)
            for (j = 0; j < 5; ++j)
                atPenaltyH[i][j] = 0.0;

        atPenaltyH[0][3] = atPenaltyH[3][0] = atp_value;
    }

    private Comparator<? super loop> loopComparator = new Comparator<loop>() {
        @Override
        public int compare(loop o, loop t1) {
            byte[] loop1 = o.loop;
            byte[] loop2 = t1.loop;
            for (int i = 0; i < loop2.length; ++i) {
                if (loop1[i] < loop2[i]) return -1;
                else if (loop1[i] > loop2[i]) return 1;
            }
            return 0;
        }
    };

    /**
     * initiates thermodynamic parameter tables of entropy and enthalpy for dimer
     */
    private void
    initMatrix() {
        int i, j;
        for (i = 1; i <= len1; ++i) {
            for (j = 1; j <= len2; ++j) {
                if (BPI[numSeq1[i]][numSeq2[j]] == 0) {
                    enthalpyDPT[i][j] = _INFINITY;
                    entropyDPT[i][j] = -1.0;
                } else {
                    enthalpyDPT[i][j] = 0.0;
                    entropyDPT[i][j] = MinEntropy;
                }
            }
        }
    }

    /**
     * initiates thermodynamic parameter tables of entropy and enthalpy for monomer
     */
    private void
    initMatrix2() {
        int i, j;
        for (i = 1; i <= len1; ++i)
            for (j = i; j <= len2; ++j) {
                if (j - i < MIN_HRPN_LOOP + 1 || (BPI[numSeq1[i]][numSeq1[j]] == 0)) {
                    enthalpyDPT[i][j] = _INFINITY;
                    entropyDPT[i][j] = -1.0;
                } else {
                    enthalpyDPT[i][j] = 0.0;
                    entropyDPT[i][j] = MinEntropy;
                }
	}
    }

    /**
     * calc-s thermod values into dynamic progr table (dimer)
     */
    private void
    fillMatrix(int maxLoop, ThalResults o) {
        int d, i, j, ii, jj;
        double[] SH;

        SH = new double[2];
        for (i = 1; i <= len1; ++i) {
            for (j = 1; j <= len2; ++j) {
                if (isFinite(enthalpyDPT[i][j])) { /* if finite */
                    SH[0] = -1.0;
                    SH[1] = _INFINITY;
                    LSH(i, j, SH);
                    if (isFinite(SH[1])) {
                        entropyDPT[i][j] = SH[0];
                        enthalpyDPT[i][j] = SH[1];
                    }
                    if (i > 1 && j > 1) {
                        maxTM(i, j); /* stack: sets entropyDPT[i][j] and enthalpyDPT[i][j] */
                        for (d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
                            ii = i - 1;
                            jj = -ii - d + (j + i);
                            if (jj < 1) {
                                ii -= abs(jj - 1);
                                jj = 1;
                            }
                            for (; ii > 0 && jj < j; --ii, ++jj) {
                                if (isFinite(enthalpyDPT[ii][jj])) {
                                    SH[0] = -1.0;
                                    SH[1] = _INFINITY;
                                    calc_bulge_internal(ii, jj, i, j, SH, 0, maxLoop);
                                    if (SH[0] < MinEntropyCutoff) {
                                        /* to not give dH any value if dS is unreasonable */
                                        SH[0] = MinEntropy;
                                        SH[1] = 0.0;
                                    }
                                    if (isFinite(SH[1])) {
                                        enthalpyDPT[i][j] = SH[1];
                                        entropyDPT[i][j] = SH[0];
                                    }
                                }
                            }
                        }
                    } /* if */
                }
            } /* for */
        } /* for */
    }

    /**
     * calc-s thermod values into dynamic progr table (monomer)
     */
    private void
    fillMatrix2(int maxLoop, ThalResults o) {
        int i, j;
        double[] SH = new double[2];
        for (j = 2; j <= len2; ++j)
            for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
                if (isFinite(enthalpyDPT[i][j])) {
                    SH[0] = -1.0;
                    SH[1] = _INFINITY;
                    maxTM2(i, j); /* calculate stack */
                    CBI(i, j, SH, 0, maxLoop); /* calculate Bulge and Internal loop and stack */
                    SH[0] = -1.0;
                    SH[1] = _INFINITY;
                    calc_hairpin(i, j, SH, 0);
                    if (isFinite(SH[1])) {
                        if (SH[0] < MinEntropyCutoff) { /* to not give dH any value if dS is unreasonable */
                            SH[0] = MinEntropy;
                            SH[1] = 0.0;
                        }
                        entropyDPT[i][j] = SH[0];
                        enthalpyDPT[i][j] = SH[1];
                    }
                }
            }
    }


    /**
     * finds max Tm while filling the dyn progr table using stacking S and stacking H (dimer)
     */
    private void
    maxTM(int i, int j) {
        double T0, T1;
        double S0, S1;
        double H0, H1;
        double[] SH = new double[2];
        T0 = T1 = -_INFINITY;
        S0 = entropyDPT[i][j];
        H0 = enthalpyDPT[i][j];
        RSH(i, j, SH);
        T0 = (H0 + dplx_init_H + SH[1]) / (S0 + dplx_init_S + SH[0] + RC); /* at current position */
        if (isFinite(enthalpyDPT[i - 1][j - 1]) && isFinite(Hs(i - 1, j - 1, 1))) {
            S1 = (entropyDPT[i - 1][j - 1] + Ss(i - 1, j - 1, 1));
            H1 = (enthalpyDPT[i - 1][j - 1] + Hs(i - 1, j - 1, 1));
            T1 = (H1 + dplx_init_H + SH[1]) / (S1 + dplx_init_S + SH[0] + RC);
        } else {
            S1 = -1.0;
            H1 = _INFINITY;
            T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
        }

        if (S1 < MinEntropyCutoff) {
            /* to not give dH any value if dS is unreasonable */
            S1 = MinEntropy;
            H1 = 0.0;
        }
        if (S0 < MinEntropyCutoff) {
            /* to not give dH any value if dS is unreasonable */
            S0 = MinEntropy;
            H0 = 0.0;
        }
        if (T1 > T0) {
            entropyDPT[i][j] = S1;
            enthalpyDPT[i][j] = H1;
        } else if (T0 >= T1) {
            entropyDPT[i][j] = S0;
            enthalpyDPT[i][j] = H0;
        }
    }

    /**
     * finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer)
     */
    private void
    maxTM2(int i, int j) {
        double T0, T1;
        double S0, S1;
        double H0, H1;
        T0 = T1 = -_INFINITY;
        S0 = entropyDPT[i][j];
        H0 = enthalpyDPT[i][j];
        T0 = (H0 + dplx_init_H) / (S0 + dplx_init_S + RC);
        if (isFinite(enthalpyDPT[i][j])) {
            S1 = (entropyDPT[i + 1][j - 1] + Ss(i, j, 2));
            H1 = (enthalpyDPT[i + 1][j - 1] + Hs(i, j, 2));
        } else {
            S1 = -1.0;
            H1 = _INFINITY;
        }
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
        if (S1 < MinEntropyCutoff) {
            S1 = MinEntropy;
            H1 = 0.0;
        }
        if (S0 < MinEntropyCutoff) {
            S0 = MinEntropy;
            H0 = 0.0;
        }

        if (T1 > T0) {
            entropyDPT[i][j] = S1;
            enthalpyDPT[i][j] = H1;
        } else {
            entropyDPT[i][j] = S0;
            enthalpyDPT[i][j] = H0;
        }
    }

    /**
     * calculate terminal entropy S and terminal enthalpy H starting reading from 5' (left) end
     */
    private void
    LSH(int i, int j, double[] EntropyEnthalpy) {
        double S1, H1, T1, G1;
        double S2, H2, T2, G2;
        S1 = S2 = -1.0;
        H1 = H2 = -_INFINITY;
        T1 = T2 = -_INFINITY;
        if (BPI[numSeq1[i]][numSeq2[j]] == 0) {
            entropyDPT[i][j] = -1.0;
            enthalpyDPT[i][j] = _INFINITY;
            return;
        }
        S1 = atPenaltyS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]][numSeq1[i - 1]];
        H1 = atPenaltyH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]][numSeq1[i - 1]];
        G1 = H1 - TEMP_KELVIN * S1;
        // TODO - Shouldn't we need to calculate T1 here?
        if (!isFinite(H1) || G1 > 0) {
            H1 = _INFINITY;
            S1 = -1.0;
            G1 = 1.0;
        }
        /** If there is two dangling ends at the same end of duplex **/
        if ((BPI[numSeq1[i - 1]][numSeq2[j - 1]] != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }
            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        } else if ((BPI[numSeq1[i - 1]][numSeq2[j - 1]] != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }
            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        } else if ((BPI[numSeq1[i - 1]][numSeq2[j - 1]] != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }
            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        }
        S2 = atPenaltyS[numSeq1[i]][numSeq2[j]];
        H2 = atPenaltyH[numSeq1[i]][numSeq2[j]];
        T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
        G1 = H1 - TEMP_KELVIN * S1;
        G2 = H2 - TEMP_KELVIN * S2;
        if (isFinite(H1)) {
            if (T1 < T2) {
                EntropyEnthalpy[0] = S2;
                EntropyEnthalpy[1] = H2;
            } else {
                EntropyEnthalpy[0] = S1;
                EntropyEnthalpy[1] = H1;
            }
        } else {
            EntropyEnthalpy[0] = S2;
            EntropyEnthalpy[1] = H2;
        }
        return;
    }

    /**
     * calculate terminal entropy S and terminal enthalpy H starting reading from 3' (right) end
     */
    private void
    RSH(int i, int j, double[] EntropyEnthalpy) {
        double G1, G2;
        double S1, S2;
        double H1, H2;
        double T1, T2;
        S1 = S2 = -1.0;
        H1 = H2 = _INFINITY;
        T1 = T2 = -_INFINITY;
        if (BPI[numSeq1[i]][numSeq2[j]] == 0) {
            EntropyEnthalpy[0] = -1.0;
            EntropyEnthalpy[1] = _INFINITY;
            return;
        }
        S1 = atPenaltyS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
        H1 = atPenaltyH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
        G1 = H1 - TEMP_KELVIN * S1;
        if (!isFinite(H1) || G1 > 0) {
            H1 = _INFINITY;
            S1 = -1.0;
            G1 = 1.0;
        }

        if (BPI[numSeq1[i + 1]][numSeq2[j + 1]] == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
                    dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
                    dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }

            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        } else if (BPI[numSeq1[i + 1]][numSeq2[j + 1]] == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }
            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        } else if (BPI[numSeq1[i + 1]][numSeq2[j + 1]] == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
            S2 = atPenaltyS[numSeq1[i]][numSeq2[j]] + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
            H2 = atPenaltyH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
            G2 = H2 - TEMP_KELVIN * S2;
            if (!isFinite(H2) || G2 > 0) {
                H2 = _INFINITY;
                S2 = -1.0;
                G2 = 1.0;
            }
            T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
            if (isFinite(H1) && G1 < 0) {
                T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
                if (T1 < T2 && G2 < 0) {
                    S1 = S2;
                    H1 = H2;
                    T1 = T2;
                }
            } else if (G2 < 0) {
                S1 = S2;
                H1 = H2;
                T1 = T2;
            }
        }
        S2 = atPenaltyS[numSeq1[i]][numSeq2[j]];
        H2 = atPenaltyH[numSeq1[i]][numSeq2[j]];
        T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
        G1 = H1 - TEMP_KELVIN * S1;
        G2 = H2 - TEMP_KELVIN * S2;
        if (isFinite(H1)) {
            if (T1 < T2) {
                EntropyEnthalpy[0] = S2;
                EntropyEnthalpy[1] = H2;
            } else {
                EntropyEnthalpy[0] = S1;
                EntropyEnthalpy[1] = H1;
            }
        } else {
            EntropyEnthalpy[0] = S2;
            EntropyEnthalpy[1] = H2;
        }
        return;
    }

    /**
     * returns stack entropy
     */
    private double Ss(int i, int j, int k) {
        if (k == 2) {
            if (i >= j)
                return -1.0;
            if (i == len1 || j == len2 + 1)
                return -1.0;

            if (i > len1)
                i -= len1;
            if (j > len2)
                j -= len2;
            return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]];
        } else {
            return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
        }
    }

    /**
     * returns stack enthalpy
     */
    private double Hs(int i, int j, int k) {
        if (k == 2) {
            if (i >= j)
                return _INFINITY;
            if (i == len1 || j == len2 + 1)
                return _INFINITY;

            if (i > len1)
                i -= len1;
            if (j > len2)
                j -= len2;
            if (isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]])) {
                return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]];
            } else {
                return _INFINITY;
            }
        } else {
            return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
        }
    }

    /**
     * carries out Bulge and Internal loop and stack calculations to hairpin
     */
    private void CBI(int i, int j, double[] EntropyEnthalpy, int traceback, int maxLoop) {
        int d, ii, jj;
        for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)
            for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
                jj = d + ii;
                if (traceback == 0) {
                    EntropyEnthalpy[0] = -1.0;
                    EntropyEnthalpy[1] = _INFINITY;
                }
                if (isFinite(enthalpyDPT[ii][jj]) && isFinite(enthalpyDPT[i][j])) {
                    calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback, maxLoop);
                    if (isFinite(EntropyEnthalpy[1])) {
                        if (EntropyEnthalpy[0] < MinEntropyCutoff) {
                            EntropyEnthalpy[0] = MinEntropy;
                            EntropyEnthalpy[1] = 0.0;
                        }
                        if (traceback == 0) {
                            enthalpyDPT[i][j] = EntropyEnthalpy[1];
                            entropyDPT[i][j] = EntropyEnthalpy[0];
                        }
                    }
                }
            }
        return;
    }

    /**
     * finds monomer structure that has maximum Tm
     */
    private void calc_hairpin(int i, int j, double[] EntropyEnthalpy, int traceback) {
        int loopSize = j - i - 1;
        double G1, G2;
        G1 = G2 = -_INFINITY;
        double[] SH = new double[2];
        SH[0] = -1.0;
        SH[1] = _INFINITY;
        if (loopSize < MIN_HRPN_LOOP) {
            EntropyEnthalpy[0] = -1.0;
            EntropyEnthalpy[1] = _INFINITY;
            return;
        }
        if (i <= len1 && len2 < j) {
            EntropyEnthalpy[0] = -1.0;
            EntropyEnthalpy[1] = _INFINITY;
            return;
        } else if (i > len2) {
            i -= len1;
            j -= len2;
        }
        if (loopSize <= 30) {
            EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
            EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
        } else {
            EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
            EntropyEnthalpy[0] = hairpinLoopEntropies[29];
        }

        if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
            EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
            EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
        } else if (loopSize == 3) { /* for loops 3 bp in length at-penalty is considered */
            EntropyEnthalpy[1] += atPenaltyH[numSeq1[i]][numSeq1[j]];
            EntropyEnthalpy[0] += atPenaltyS[numSeq1[i]][numSeq1[j]];
        }

        if (loopSize == 3) {         /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
            if (numTriloops > 0) {
                loop key = new loop(Arrays.copyOfRange(numSeq1, i, i + 5));
                int index = -1;
                if ((index = binarySearch(triloopEnthalpies, key, loopComparator)) >= 0)
                    EntropyEnthalpy[1] += triloopEnthalpies[index].value;
                if ((index = binarySearch(triloopEntropies, key, loopComparator)) >= 0)
                    EntropyEnthalpy[0] += triloopEntropies[index].value;
            }
        } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
            if (numTetraloops > 0) {
                loop key = new loop(Arrays.copyOfRange(numSeq1, i, i + 6));
                int index = -1;
                if ((index = binarySearch(tetraloopEnthalpies, key, loopComparator)) >= 0)
                    EntropyEnthalpy[1] += tetraloopEnthalpies[index].value;
                if ((index = binarySearch(tetraloopEntropies, key, loopComparator)) >= 0)
                    EntropyEnthalpy[0] += tetraloopEntropies[index].value;
            }
        }
        if (!isFinite(EntropyEnthalpy[1])) {
            EntropyEnthalpy[1] = _INFINITY;
            EntropyEnthalpy[0] = -1.0;
        }
        if (isPositive(EntropyEnthalpy[1]) && isPositive(EntropyEnthalpy[0]) && (!isPositive(enthalpyDPT[i][j]) || !isPositive(entropyDPT[i][j]))) { /* if both, S and H are positive */
            EntropyEnthalpy[1] = _INFINITY;
            EntropyEnthalpy[0] = -1.0;
        }
        RSH(i, j, SH);
        G1 = EntropyEnthalpy[1] + SH[1] - TEMP_KELVIN * (EntropyEnthalpy[0] + SH[0]);
        G2 = enthalpyDPT[i][j] + SH[1] - TEMP_KELVIN * (entropyDPT[i][j] + SH[0]);
        if (G2 < G1 && traceback == 0) {
            EntropyEnthalpy[0] = entropyDPT[i][j];
            EntropyEnthalpy[1] = enthalpyDPT[i][j];
        }
        return;
    }


    /**
     * calculates bulges and internal loops for dimer structures
     */
    private void
    calc_bulge_internal(int i, int j, int ii, int jj, double[] EntropyEnthalpy, int traceback, int maxLoop) {
        int loopSize1, loopSize2, loopSize;
        double S, H, G1, G2;
        int N, N_loop;
        double[] SH = new double[2];
        SH[0] = -1.0;
        SH[1] = _INFINITY;
        S = -1.0;
        H = _INFINITY;
        loopSize1 = ii - i - 1;
        loopSize2 = jj - j - 1;
        if (ii < jj) {
            N = ((2 * i) / 2);
            N_loop = N;
            if (loopSize1 > 2) N_loop -= (loopSize1 - 2);
            if (loopSize2 > 2) N_loop -= (loopSize2 - 2);
        } else {
            N = ((2 * j) / 2);
            N_loop = 2 * jj;
            if (loopSize1 > 2) N_loop -= (loopSize1 - 2);
            if (loopSize2 > 2) N_loop -= (loopSize2 - 2);
            N_loop = (N_loop / 2) - 1;
        }
        if (DEBUG) {
            if (ii <= i) {
                System.err.print("Error in calc_bulge_internal(): ii is not greater than i\n");
            }
            if (jj <= j)
                System.err.print("Error in calc_bulge_internal(): jj is not greater than j\n");

            if (loopSize1 + loopSize2 > maxLoop) {
                System.err.print("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n");
                return;
            }
            if (loopSize1 == 0 && loopSize2 == 0) {
                System.err.print("Error: calc_bulge_internal() called with nonsense\n");
                return;
            }
        }
        loopSize = loopSize1 + loopSize2 - 1;
        if ((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
            if (loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                                  the intervening nn-pair must be added */

                if ((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
                    H = bulgeLoopEnthalpies[loopSize] +
                            stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
                    S = bulgeLoopEntropies[loopSize] +
                            stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
                }
                if (isPositive(H) || isPositive(S)) {
                    H = _INFINITY;
                    S = -1.0;
                }
                H += enthalpyDPT[i][j];
                S += entropyDPT[i][j];
                if (!isFinite(H)) {
                    H = _INFINITY;
                    S = -1.0;
                }
                RSH(ii, jj, SH);
                G1 = H + SH[1] - TEMP_KELVIN * (S + SH[0]);
                G2 = enthalpyDPT[ii][jj] + SH[1] - TEMP_KELVIN * ((entropyDPT[ii][jj] + SH[0]));
                if ((G1 < G2) || (traceback == 1)) {
                    EntropyEnthalpy[0] = S;
                    EntropyEnthalpy[1] = H;
                }
            } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

                H = bulgeLoopEnthalpies[loopSize] + atPenaltyH[numSeq1[i]][numSeq2[j]] + atPenaltyH[numSeq1[ii]][numSeq2[jj]];
                H += enthalpyDPT[i][j];

                S = bulgeLoopEntropies[loopSize] + atPenaltyS[numSeq1[i]][numSeq2[j]] + atPenaltyS[numSeq1[ii]][numSeq2[jj]];
                S += entropyDPT[i][j];
                if (!isFinite(H)) {
                    H = _INFINITY;
                    S = -1.0;
                }
                if (isPositive(H) && isPositive(S)) {
                    H = _INFINITY;
                    S = -1.0;
                }

                RSH(ii, jj, SH);
                G1 = H + SH[1] - TEMP_KELVIN * (S + SH[0]);
                G2 = enthalpyDPT[ii][jj] + SH[1] - TEMP_KELVIN * (entropyDPT[ii][jj] + SH[0]);
                if (G1 < G2 || (traceback == 1)) {
                    EntropyEnthalpy[0] = S;
                    EntropyEnthalpy[1] = H;
                }

            }
        } else if (loopSize1 == 1 && loopSize2 == 1) {
            S = stackint2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]] +
                    stackint2Entropies[numSeq2[jj]][numSeq2[jj - 1]][numSeq1[ii]][numSeq1[ii - 1]];
            S += entropyDPT[i][j];

            H = stackint2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]] +
                    stackint2Enthalpies[numSeq2[jj]][numSeq2[jj - 1]][numSeq1[ii]][numSeq1[ii - 1]];
            H += enthalpyDPT[i][j];
            if (!isFinite(H)) {
                H = _INFINITY;
                S = -1.0;
            }
            if (isPositive(H) && isPositive(S)) {
                H = _INFINITY;
                S = -1.0;
            }
            RSH(ii, jj, SH);
            G1 = H + SH[1] - TEMP_KELVIN * (S + SH[0]);
            G2 = enthalpyDPT[ii][jj] + SH[1] - TEMP_KELVIN * (entropyDPT[ii][jj] + SH[0]);
            if ((G1 < G2) || traceback == 1) {
                EntropyEnthalpy[0] = S;
                EntropyEnthalpy[1] = H;
            }
            return;
        } else { /* only internal loops */
            H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]] +
                    tstackEnthalpies[numSeq2[jj]][numSeq2[jj - 1]][numSeq1[ii]][numSeq1[ii - 1]]
                    + (ILAH * abs(loopSize1 - loopSize2));
            H += enthalpyDPT[i][j];

            S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]] +
                    tstackEntropies[numSeq2[jj]][numSeq2[jj - 1]][numSeq1[ii]][numSeq1[ii - 1]] + (ILAS * abs(loopSize1 - loopSize2));
            S += entropyDPT[i][j];
            if (!isFinite(H)) {
                H = _INFINITY;
                S = -1.0;
            }
            if (isPositive(H) && isPositive(S)) {
                H = _INFINITY;
                S = -1.0;
            }
            RSH(ii, jj, SH);
            G1 = H + SH[1] - TEMP_KELVIN * (S + SH[0]);
            G2 = enthalpyDPT[ii][jj] + SH[1] - TEMP_KELVIN * (entropyDPT[ii][jj] + SH[0]);
            if ((G1 < G2) || (traceback == 1)) {
                EntropyEnthalpy[0] = S;
                EntropyEnthalpy[1] = H;
            }
        }
        return;
    }

    /**
     * calculates bulges and internal loops for monomer structures
     */
    private void
    calc_bulge_internal2(int i, int j, int ii, int jj, double[] EntropyEnthalpy, int traceback, int maxLoop) {
        int loopSize1, loopSize2, loopSize;
        double T1, T2;
        double S, H;
        /* int N, N_loop; Triinu, please review */
        T1 = T2 = -_INFINITY;
        S = MinEntropy;
        H = 0.0;
        loopSize1 = ii - i - 1;
        loopSize2 = j - jj - 1;
        if (loopSize1 + loopSize2 > maxLoop) {
            EntropyEnthalpy[0] = -1.0;
            EntropyEnthalpy[1] = _INFINITY;
            return;
        }
        /* Triinu, please review the statements below. */
        /* if(i < (len1 -j)) { */
        /* N  = i; */
        /* N_loop = (i - 1); */
        /* } else { */
        /* N = len1-j;  */
        /* N_loop = len1 - j - 1; */
        /* } */
        if (DEBUG) {
            if (ii <= i)
                System.err.print("Error in calc_bulge_internal(): ii isn't greater than i\n");
            if (jj >= j)
                System.err.print("Error in calc_bulge_internal(): jj isn't less than j\n");
            if (ii >= jj)
                System.err.print("Error in calc_bulge_internal(): jj isn't greater than ii\n");

            if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j)) {
                EntropyEnthalpy[0] = -1.0;
                EntropyEnthalpy[1] = _INFINITY;
                return;
            }

            if (loopSize1 + loopSize2 > maxLoop) {
                System.err.print("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n");
                return;
            }
            if (loopSize1 == 0 && loopSize2 == 0) {
                System.err.print("Error: calc_bulge_internal() called with nonsense\n");
                return;
            }

            if (i > len1)
                i -= len1;
            if (ii > len1)
                ii -= len1;
            if (j > len2)
                j -= len2;
            if (jj > len2)
                jj -= len2;
        }
        loopSize = loopSize1 + loopSize2 - 1; /* for indx only */
        if ((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
            if (loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                                  the intervening nn-pair must be added */
                if ((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
                    H = bulgeLoopEnthalpies[loopSize] +
                            stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
                    S = bulgeLoopEntropies[loopSize] +
                            stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
                }
                if (traceback != 1) {
                    H += enthalpyDPT[ii][jj]; /* bulge koos otsaga, st bulge i,j-ni */
                    S += entropyDPT[ii][jj];
                }

                if (!isFinite(H)) {
                    H = _INFINITY;
                    S = -1.0;
                }

                T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
                T2 = (enthalpyDPT[i][j] + dplx_init_H) / ((entropyDPT[i][j]) + dplx_init_S + RC);

                if ((T1 > T2) || ((traceback != 0 && T1 >= T2) || traceback == 1)) {
                    EntropyEnthalpy[0] = S;
                    EntropyEnthalpy[1] = H;
                }

            } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

                H = bulgeLoopEnthalpies[loopSize] + atPenaltyH[numSeq1[i]][numSeq2[j]] + atPenaltyH[numSeq1[ii]][numSeq2[jj]];
                if (traceback != 1)
                    H += enthalpyDPT[ii][jj];

                S = bulgeLoopEntropies[loopSize] + atPenaltyS[numSeq1[i]][numSeq2[j]] + atPenaltyS[numSeq1[ii]][numSeq2[jj]];
                if (traceback != 1)
                    S += entropyDPT[ii][jj];
                if (!isFinite(H)) {
                    H = _INFINITY;
                    S = -1.0;
                }

                T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
                T2 = (enthalpyDPT[i][j] + dplx_init_H) / (entropyDPT[i][j] + dplx_init_S + RC);

                if ((T1 > T2) || ((traceback != 0 && T1 >= T2) || (traceback == 1))) {
                    EntropyEnthalpy[0] = S;
                    EntropyEnthalpy[1] = H;
                }
            }
        } /* end of calculating bulges */ else if (loopSize1 == 1 && loopSize2 == 1) {
            /* mismatch nearest neighbor parameters */

            S = stackint2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]] +
                    stackint2Entropies[numSeq2[jj]][numSeq2[jj + 1]][numSeq1[ii]][numSeq1[ii - 1]];
            if (traceback != 1)
                S += entropyDPT[ii][jj];

            H = stackint2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]] +
                    stackint2Enthalpies[numSeq2[jj]][numSeq2[jj + 1]][numSeq1[ii]][numSeq1[ii - 1]];
            if (traceback != 1)
                H += enthalpyDPT[ii][jj];
            if (!isFinite(H)) {
                H = _INFINITY;
                S = -1.0;
            }

            T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
            T2 = (enthalpyDPT[i][j] + dplx_init_H) / (entropyDPT[i][j] + dplx_init_S + RC);
            if ((DBL_EQ(T1, T2) == 2) || traceback != 0) {
                if ((T1 > T2) || ((traceback != 0 && T1 >= T2) || traceback == 1)) {
                    EntropyEnthalpy[0] = S;
                    EntropyEnthalpy[1] = H;
                }
            }
            return;
        } else { /* only internal loops */

            H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]] +
                    tstackEnthalpies[numSeq2[jj]][numSeq2[jj + 1]][numSeq1[ii]][numSeq1[ii - 1]]
                    + (ILAH * abs(loopSize1 - loopSize2));
            if (traceback != 1)
                H += enthalpyDPT[ii][jj];

            S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]] +
                    tstackEntropies[numSeq2[jj]][numSeq2[jj + 1]][numSeq1[ii]][numSeq1[ii - 1]] + (ILAS * abs(loopSize1 - loopSize2));
            if (traceback != 1)
                S += entropyDPT[ii][jj];
            if (!isFinite(H)) {
                H = _INFINITY;
                S = -1.0;
            }

            T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
            T2 = (enthalpyDPT[i][j] + dplx_init_H) / ((entropyDPT[i][j]) + dplx_init_S + RC);
            if ((T1 > T2) || ((traceback != 0 && T1 >= T2) || (traceback == 1))) {
                EntropyEnthalpy[0] = S;
                EntropyEnthalpy[1] = H;
            }
        }
        return;
    }

    /**
     * terminal bp for monomer structure
     */
    private void
    calc_terminal_bp(double temp) { /* compute exterior loop */
        int i;
        int max;
        send5[0] = send5[1] = -1.0;
        hend5[0] = hend5[1] = _INFINITY;
        for (i = 2; i <= (len1); i++) {
            send5[i] = MinEntropy;
            hend5[i] = 0;
        }

        double T1, T2, T3, T4, T5;
        T1 = T2 = T3 = T4 = T5 = -_INFINITY;
        double G;
        /* adding terminal penalties to 3' end and to 5' end */
        for (i = 2; i <= len1; ++i) {
            max = 0;
            T1 = T2 = T3 = T4 = T5 = -_INFINITY;
            T1 = (hend5[i - 1] + dplx_init_H) / (send5[i - 1] + dplx_init_S + RC);
            T2 = (END5_1(i, 1) + dplx_init_H) / (END5_1(i, 2) + dplx_init_S + RC);
            T3 = (END5_2(i, 1) + dplx_init_H) / (END5_2(i, 2) + dplx_init_S + RC);
            T4 = (END5_3(i, 1) + dplx_init_H) / (END5_3(i, 2) + dplx_init_S + RC);
            T5 = (END5_4(i, 1) + dplx_init_H) / (END5_4(i, 2) + dplx_init_S + RC);
            max = max5(T1, T2, T3, T4, T5);
            switch (max) {
                case 1:
                    send5[i] = send5[i - 1];
                    hend5[i] = hend5[i - 1];
                    break;
                case 2:
                    G = END5_1(i, 1) - (temp * (END5_1(i, 2)));
                    if (G < G2) {
                        send5[i] = END5_1(i, 2);
                        hend5[i] = END5_1(i, 1);
                    } else {
                        send5[i] = send5[i - 1];
                        hend5[i] = hend5[i - 1];
                    }
                    break;
                case 3:
                    G = END5_2(i, 1) - (temp * (END5_2(i, 2)));
                    if (G < G2) {
                        send5[i] = END5_2(i, 2);
                        hend5[i] = END5_2(i, 1);
                    } else {
                        send5[i] = send5[i - 1];
                        hend5[i] = hend5[i - 1];
                    }
                    break;
                case 4:
                    G = END5_3(i, 1) - (temp * (END5_3(i, 2)));
                    if (G < G2) {
                        send5[i] = END5_3(i, 2);
                        hend5[i] = END5_3(i, 1);
                    } else {
                        send5[i] = send5[i - 1];
                        hend5[i] = hend5[i - 1];
                    }
                    break;
                case 5:
                    G = END5_4(i, 1) - (temp * (END5_4(i, 2)));
                    if (G < G2) {
                        send5[i] = END5_4(i, 2);
                        hend5[i] = END5_4(i, 1);
                    } else {
                        send5[i] = send5[i - 1];
                        hend5[i] = hend5[i - 1];
                    }
                    break;
                default:
                    if (DEBUG) System.err.printf("WARNING: max5 returned character code %d ??\n", max);
                    break;
            }
        }
    }


    /** executed in calc_terminal_bp; to find structure that corresponds to max Tm for terminal bp */
    /**
     * END5_1(X,1/2) - 1=Enthalpy, 2=Entropy
     */
    private double
    END5_1(int i, int hs) {
        int k;
        double max_tm; /* energy min */
        double T1, T2;
        double H, S;
        double H_max, S_max;
        H_max = H = _INFINITY;
        S_max = S = -1.0;
        T1 = T2 = -_INFINITY;
        max_tm = -_INFINITY;
        for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
            T1 = (hend5[k] + dplx_init_H) / (send5[k] + dplx_init_S + RC);
            T2 = (0 + dplx_init_H) / (0 + dplx_init_S + RC);
            if (T1 >= T2) {
                H = hend5[k] + atPenaltyH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i];
                S = send5[k] + atPenaltyS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i];
                if (!isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            } else {
                H = 0 + atPenaltyH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i];
                S = 0 + atPenaltyS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            }
            if (max_tm < T1) {
                if (S > MinEntropyCutoff) {
                    H_max = H;
                    S_max = S;
                    max_tm = T1;
                }
            }
        }
        if (hs == 1) return H_max;
        return S_max;
    }

    private double
    END5_2(int i, int hs) {
        int k;
        double max_tm;
        double T1, T2;
        double H, S;
        double H_max, S_max;
        H_max = H = _INFINITY;
        T1 = T2 = max_tm = -_INFINITY;
        S_max = S = -1.0;
        for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
            T1 = (hend5[k] + dplx_init_H) / (send5[k] + dplx_init_S + RC);
            T2 = (0 + dplx_init_H) / (0 + dplx_init_S + RC);
            if (T1 >= T2) {
                H = hend5[k] + atPenaltyH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + enthalpyDPT[k + 2][i];
                S = send5[k] + atPenaltyS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + entropyDPT[k + 2][i];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            } else {
                H = 0 + atPenaltyH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + enthalpyDPT[k + 2][i];
                S = 0 + atPenaltyS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + entropyDPT[k + 2][i];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            }
            if (max_tm < T1) {
                if (S > MinEntropyCutoff) {
                    H_max = H;
                    S_max = S;
                    max_tm = T1;
                }
            }
        }
        if (hs == 1) return H_max;
        return S_max;
    }

    private double
    END5_3(int i, int hs) {
        int k;
        double max_tm;
        double T1, T2;
        double H, S;
        double H_max, S_max;
        H_max = H = _INFINITY;
        ;
        T1 = T2 = max_tm = -_INFINITY;
        S_max = S = -1.0;
        for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
            T1 = (hend5[k] + dplx_init_H) / (send5[k] + dplx_init_S + RC);
            T2 = (0 + dplx_init_H) / (0 + dplx_init_S + RC);
            if (T1 >= T2) {
                H = hend5[k] + atPenaltyH[numSeq1[k + 1]][numSeq1[i - 1]] + Hd3(i - 1, k + 1) + enthalpyDPT[k + 1][i - 1];
                S = send5[k] + atPenaltyS[numSeq1[k + 1]][numSeq1[i - 1]] + Sd3(i - 1, k + 1) + entropyDPT[k + 1][i - 1];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            } else {
                H = 0 + atPenaltyH[numSeq1[k + 1]][numSeq1[i - 1]] + Hd3(i - 1, k + 1) + enthalpyDPT[k + 1][i - 1];
                S = 0 + atPenaltyS[numSeq1[k + 1]][numSeq1[i - 1]] + Sd3(i - 1, k + 1) + entropyDPT[k + 1][i - 1];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            }
            if (max_tm < T1) {
                if (S > MinEntropyCutoff) {
                    H_max = H;
                    S_max = S;
                    max_tm = T1;
                }
            }
        }
        if (hs == 1) return H_max;
        return S_max;
    }

    private double
    END5_4(int i, int hs) {
        int k;
        double max_tm;
        double T1, T2;
        double H, S;
        double H_max, S_max;
        H_max = H = _INFINITY;
        T1 = T2 = max_tm = -_INFINITY;
        S_max = S = -1.0;
        for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
            T1 = (hend5[k] + dplx_init_H) / (send5[k] + dplx_init_S + RC);
            T2 = (0 + dplx_init_H) / (0 + dplx_init_S + RC);
            if (T1 >= T2) {
                H = hend5[k] + atPenaltyH[numSeq1[k + 2]][numSeq1[i - 1]] + Htstack(i - 1, k + 2) + enthalpyDPT[k + 2][i - 1];
                S = send5[k] + atPenaltyS[numSeq1[k + 2]][numSeq1[i - 1]] + Ststack(i - 1, k + 2) + entropyDPT[k + 2][i - 1];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            } else {
                H = 0 + atPenaltyH[numSeq1[k + 2]][numSeq1[i - 1]] + Htstack(i - 1, k + 2) + enthalpyDPT[k + 2][i - 1];
                S = 0 + atPenaltyS[numSeq1[k + 2]][numSeq1[i - 1]] + Ststack(i - 1, k + 2) + entropyDPT[k + 2][i - 1];
                if (!isFinite(H) || H > 0 || S > 0) {
                    H = _INFINITY;
                    S = -1.0;
                }
                T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
            }
            if (max_tm < T1) {
                if (S > MinEntropyCutoff) {
                    H_max = H;
                    S_max = S;
                    max_tm = T1;
                }
            }
        }
        if (hs == 1) return H_max;
        return S_max;
    }


    /**
     * returns thermodynamic value (S) for 5' dangling end
     */
    private double Sd5(int i, int j) {
        return dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
    }

    /**
     * returns thermodynamic value (H) for 5' dangling end
     */
    private double Hd5(int i, int j) {
        return dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
    }

    /**
     * returns thermodynamic value (S) for 3' dangling end
     */
    private double Sd3(int i, int j) {
        return dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]];
    }

    /**
     * returns thermodynamic value (H) for 3' dangling end
     */
    private double Hd3(int i, int j) {
        return dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]];
    }

    /**
     * returns entropy value for terminal stack
     */
    private double Ststack(int i, int j) {
        return tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
    }

    /**
     * returns enthalpy value for terminal stack
     */
    private double Htstack(int i, int j) { /* e.g AG_TC 210 */
        return tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
    }

    /**
     * Return true if string is symmetrical, false otherwise.
     */
    private boolean symmetry_thermo(String seq_str) {
        char s;
        char e;
        int i = 0;
        int seq_len = seq_str.length();
        int seq = 0, seq_end = 0;
        ;
        int mp = seq_len / 2;
        if (seq_len % 2 == 1) {
            return false;
        }
        seq_end += seq_len;
        seq_end--;
        seq_str = seq_str.toUpperCase();
        while (i < mp) {
            i++;
            s = seq_str.charAt(seq);
            e = seq_str.charAt(seq_end);
            if ((s == 'A' && e != 'T')
                    || (s == 'T' && e != 'A')
                    || (e == 'A' && s != 'T')
                    || (e == 'T' && s != 'A')) {
                return false;
            }
            if ((s == 'C' && e != 'G')
                    || (s == 'G' && e != 'C')
                    || (e == 'C' && s != 'G')
                    || (e == 'G' && s != 'C')) {
                return false;
            }
            seq++;
            seq_end--;
        }
        return true;
    }

    /**
     * returns length of char; to avoid warnings while compiling
     */
    private int length_unsig_char(String str) {
        return str.length();
    }

    /**
     * traceback for  unimolecular structure (hairpins)
     */
    private void tracebacku(int[] bp, int maxLoop, ThalResults o) {
        int i, j;
        i = j = 0;
        int ii, jj, k;
        tracer top, stack = null;
        double[] SH1 = new double[2];
        double[] SH2 = new double[2];
        double[] EntropyEnthalpy = new double[2];
        stack = push(stack, len1, 0, 1);
        while (stack != null) {
            top = stack;
            stack = stack.next;
            i = top.i;
            j = top.j;
            if (top.mtrx == 1) {
                while (equal(send5[i], send5[i - 1]) && equal(hend5[i], hend5[i - 1])) /* if previous structure is the same as this one */
                    --i;
                if (i == 0)
                    continue;
                if (equal(send5[i], END5_1(i, 2)) && equal(hend5[i], END5_1(i, 1))) {
                    for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
                        if (equal(send5[i], atPenaltyS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i]) &&
                                equal(hend5[i], atPenaltyH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i])) {
                            stack = push(stack, k + 1, i, 0);
                            break;
                        } else if (equal(send5[i], send5[k] + atPenaltyS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i]) &&
                                equal(hend5[i], hend5[k] + atPenaltyH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i])) {
                            stack = push(stack, k + 1, i, 0);
                            stack = push(stack, k, 0, 1);
                            break;
                        }
                } else if (equal(send5[i], END5_2(i, 2)) && equal(hend5[i], END5_2(i, 1))) {
                    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
                        if (equal(send5[i], atPenaltyS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + entropyDPT[k + 2][i]) &&
                                equal(hend5[i], atPenaltyH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + enthalpyDPT[k + 2][i])) {
                            stack = push(stack, k + 2, i, 0);
                            break;
                        } else if (equal(send5[i], send5[k] + atPenaltyS[numSeq1[k + 2]][numSeq1[i]] + Sd5(i, k + 2) + entropyDPT[k + 2][i]) &&
                                equal(hend5[i], hend5[k] + atPenaltyH[numSeq1[k + 2]][numSeq1[i]] + Hd5(i, k + 2) + enthalpyDPT[k + 2][i])) {
                            stack = push(stack, k + 2, i, 0);
                            stack = push(stack, k, 0, 1);
                            break;
                        }
                } else if (equal(send5[i], END5_3(i, 2)) && equal(hend5[i], END5_3(i, 1))) {
                    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
                        if (equal(send5[i], atPenaltyS[numSeq1[k + 1]][numSeq1[i - 1]] + Sd3(i - 1, k + 1) + entropyDPT[k + 1][i - 1])
                                && equal(hend5[i], atPenaltyH[numSeq1[k + 1]][numSeq1[i - 1]] + Hd3(i - 1, k + 1) + enthalpyDPT[k + 1][i - 1])) {
                            stack = push(stack, k + 1, i - 1, 0);
                            break;
                        } else if (equal(send5[i], send5[k] + atPenaltyS[numSeq1[k + 1]][numSeq1[i - 1]] + Sd3(i - 1, k + 1) + entropyDPT[k + 1][i - 1]) &&
                                equal(hend5[i], hend5[k] + atPenaltyH[numSeq1[k + 1]][numSeq1[i - 1]] + Hd3(i - 1, k + 1) + enthalpyDPT[k + 1][i - 1])) {
                            stack = push(stack, k + 1, i - 1, 0); /* matrix 0  */
                            stack = push(stack, k, 0, 1); /* matrix 3 */
                            break;
                        }
                } else if (equal(send5[i], END5_4(i, 2)) && equal(hend5[i], END5_4(i, 1))) {
                    for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
                        if (equal(send5[i], atPenaltyS[numSeq1[k + 2]][numSeq1[i - 1]] + Ststack(i - 1, k + 2) + entropyDPT[k + 2][i - 1]) &&
                                equal(hend5[i], atPenaltyH[numSeq1[k + 2]][numSeq1[i - 1]] + Htstack(i - 1, k + 2) + enthalpyDPT[k + 2][i - 1])) {
                            stack = push(stack, k + 2, i - 1, 0);
                            break;
                        } else if (equal(send5[i], send5[k] + atPenaltyS[numSeq1[k + 2]][numSeq1[i - 1]] + Ststack(i - 1, k + 2) + entropyDPT[k + 2][i - 1]) &&
                                equal(hend5[i], hend5[k] + atPenaltyH[numSeq1[k + 2]][numSeq1[i - 1]] + Htstack(i - 1, k + 2) + enthalpyDPT[k + 2][i - 1])) {
                            stack = push(stack, k + 2, i - 1, 0);
                            stack = push(stack, k, 0, 1);
                            break;
                        }
                }
            } else if (top.mtrx == 0) {
                bp[i - 1] = j;
                bp[j - 1] = i;
                SH1[0] = -1.0;
                SH1[1] = _INFINITY;
                calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
                SH2[0] = -1.0;
                SH2[1] = _INFINITY;
                CBI(i, j, SH2, 2, maxLoop);
                if (equal(entropyDPT[i][j], Ss(i, j, 2) + entropyDPT[i + 1][j - 1]) &&
                        equal(enthalpyDPT[i][j], Hs(i, j, 2) + enthalpyDPT[i + 1][j - 1])) {
                    stack = push(stack, i + 1, j - 1, 0);
                } else if (equal(entropyDPT[i][j], SH1[0]) && equal(enthalpyDPT[i][j], SH1[1])) ;
                else if (equal(entropyDPT[i][j], SH2[0]) && equal(enthalpyDPT[i][j], SH2[1])) {
                    int d, done;
                    for (done = 0, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && done == 0; --d)
                        for (ii = i + 1; ii < j - d; ++ii) {
                            jj = d + ii;
                            EntropyEnthalpy[0] = -1.0;
                            EntropyEnthalpy[1] = _INFINITY;
                            calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, 1, maxLoop);
                            if (equal(entropyDPT[i][j], EntropyEnthalpy[0] + entropyDPT[ii][jj]) &&
                                    equal(enthalpyDPT[i][j], EntropyEnthalpy[1] + enthalpyDPT[ii][jj])) {
                                stack = push(stack, ii, jj, 0);
                                ++done;
                                break;
                            }
                        }
                } else {
                }
            }
        }
    }

    /**
     * traceback for dimers
     */
    private void traceback(int i, int j, double RT, int[] ps1, int[] ps2, int maxLoop, ThalResults o) {
        int d, ii, jj;
        boolean done;
        double[] SH = new double[2];
        ps1[i - 1] = j;
        ps2[j - 1] = i;
        while (true) {
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            LSH(i, j, SH);
            if (equal(entropyDPT[i][j], SH[0]) && equal(enthalpyDPT[i][j], SH[1])) {
                break;
            }
            done = false;
            if (i > 1 && j > 1 && equal(entropyDPT[i][j], Ss(i - 1, j - 1, 1) + entropyDPT[i - 1][j - 1]) && equal(enthalpyDPT[i][j], Hs(i - 1, j - 1, 1) + enthalpyDPT[i - 1][j - 1])) {
                i = i - 1;
                j = j - 1;
                ps1[i - 1] = j;
                ps2[j - 1] = i;
                done = true;
            }
            for (d = 3; !done && d <= maxLoop + 2; ++d) {
                ii = i - 1;
                jj = -ii - d + (j + i);
                if (jj < 1) {
                    ii -= abs(jj - 1);
                    jj = 1;
                }
                for (; !done && ii > 0 && jj < j; --ii, ++jj) {
                    SH[0] = -1.0;
                    SH[1] = _INFINITY;
                    calc_bulge_internal(ii, jj, i, j, SH, 1, maxLoop);
                    if (equal(entropyDPT[i][j], SH[0]) && equal(enthalpyDPT[i][j], SH[1])) {
                        i = ii;
                        j = jj;
                        ps1[i - 1] = j;
                        ps2[j - 1] = i;
                        done = true;
                        break;
                    }
                }
            }
        }
    }

    /**
     * prints ascii output of dimer structure
     */
    String
    drawDimer(int[] ps1, int[] ps2, double temp, double H, double S, final int mode, double t37, ThalResults o) {
        int ret_nr, ret_pr_once;
        String ret_para = "";
        int i, j, k, numSS1, numSS2, N;
        double G, t;
        t = G = 0;
        if (!isFinite(temp)) {
            o.temp = 0.0;
            o.msg = "Could not compute thermodynamic values";
            return "No predicted secondary structures for given sequences";
        } else {
            N = 0;
            for (i = 0; i < len1; i++) if (ps1[i] > 0) ++N;
            for (i = 0; i < len2; i++) if (ps2[i] > 0) ++N;
            N = (N / 2) - 1;
            S = S + (N * saltCorrection);
            t = (H / (S + RC)) - ABSOLUTE_ZERO;
            G = H - (t37 * S);
            o.temp = t;
            o.delta_G = G;
            o.delta_H = H;
            o.delta_S = S;
            o.msg = String.format("Tm: %,12.2f \u00B0C\n\u0394G: %,12.2f cal/mol\n\u0394H: %,12.2f cal/mol\n\u0394S: %,12.2f cal/mol*K\n", t, G, H, S);
            if (mode == OUTPUT_MODE_FAST) return null;
        }

        StringBuilder[] duplex = new StringBuilder[4];
        for (int ii = 0; ii < 4; ii++) duplex[ii] = new StringBuilder();

        i = 0;
        numSS1 = 0;
        while (ps1[i++] == 0) ++numSS1;
        j = 0;
        numSS2 = 0;
        while (ps2[j++] == 0) ++numSS2;

        if (numSS1 >= numSS2) {
            for (i = 0; i < numSS1; ++i) {
                append(duplex[0], oligo1.charAt(i));
                append(duplex[1], ' ');
                append(duplex[2], ' ');
            }
            for (j = 0; j < numSS1 - numSS2; ++j) append(duplex[3], ' ');
            for (j = 0; j < numSS2; ++j) append(duplex[3], oligo2.charAt(j));
        } else {
            for (j = 0; j < numSS2; ++j) {
                append(duplex[3], oligo2.charAt(j));
                append(duplex[1], ' ');
                append(duplex[2], ' ');
            }
            for (i = 0; i < numSS2 - numSS1; ++i)
                append(duplex[0], ' ');
            for (i = 0; i < numSS1; ++i)
                append(duplex[0], oligo1.charAt(i));
        }
        i = numSS1 + 1;
        j = numSS2 + 1;

        while (i <= len1) {
            while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
                append(duplex[0], ' ');
                append(duplex[1], oligo1.charAt(i - 1));
                append(duplex[2], oligo2.charAt(j - 1));
                append(duplex[3], ' ');
                ++i;
                ++j;
            }
            numSS1 = 0;
            while (i <= len1 && ps1[i - 1] == 0) {
                append(duplex[0], oligo1.charAt(i - 1));
                append(duplex[1], ' ');
                ++numSS1;
                ++i;
            }
            numSS2 = 0;
            while (j <= len2 && ps2[j - 1] == 0) {
                append(duplex[2], ' ');
                append(duplex[3], oligo2.charAt(j - 1));
                ++numSS2;
                ++j;
            }
            if (numSS1 < numSS2)
                for (k = 0; k < numSS2 - numSS1; ++k) {
                    append(duplex[0], '-');
                    append(duplex[1], ' ');
                }
            else if (numSS1 > numSS2)
                for (k = 0; k < numSS1 - numSS2; ++k) {
                    append(duplex[2], ' ');
                    append(duplex[3], '-');
                }
        }
        if ((mode == OUTPUT_MODE_TEXT)) {
            return String.format("SEQ\t%s\nSEQ\t%s\nSTR\t%s\nSTR\t%s\n", duplex[0], duplex[1], duplex[2], duplex[3]);
        }
        if (mode == OUTPUT_MODE_GRAPHICAL) {
            StringBuilder[] ret_str = new StringBuilder[4];
            for (int ii = 0; ii < 4; ii++) ret_str[ii] = new StringBuilder();
            /* Join top primer */
            strcpy(ret_str[0], "   ");
            strcat(ret_str[0], duplex[0]);
            ret_nr = 0;
            while (ret_nr < duplex[1].length()) {
                if (duplex[1].charAt(ret_nr) == 'A' || duplex[1].charAt(ret_nr) == 'T' ||
                        duplex[1].charAt(ret_nr) == 'C' || duplex[1].charAt(ret_nr) == 'G' ||
                        duplex[1].charAt(ret_nr) == '-') {
                    replace(ret_str[0], ret_nr + 3, duplex[1].charAt(ret_nr));
                }
                ret_nr++;
            }
            if (strlen(duplex[1]) > strlen(duplex[0])) {
                truncate(ret_str[0], strlen(duplex[1]) + 3);
            }
            /* Clean Ends */
            ret_nr = strlen(ret_str[0]) - 1;
            while (ret_nr > 0 && (ret_str[0].charAt(ret_nr) == ' ' || ret_str[0].charAt(ret_nr) == '-')) {
                truncate(ret_str[0], ret_nr--);
            }
            /* Write the 5' */
            ret_nr = 3;
            ret_pr_once = 1;
            while (ret_nr < ret_str[0].length() && ret_pr_once == 1) {
                if (ret_str[0].charAt(ret_nr) == 'A' || ret_str[0].charAt(ret_nr) == 'T' ||
                        ret_str[0].charAt(ret_nr) == 'C' || ret_str[0].charAt(ret_nr) == 'G' ||
                        ret_str[0].charAt(ret_nr) == '-') {
                    replace(ret_str[0], ret_nr - 3, '5');
                    replace(ret_str[0], ret_nr - 2, '\'');
                    ret_pr_once = 0;
                }
                ret_nr++;
            }

            /* Create the align tics */
            strcpy(ret_str[1], "     ");
            for (i = 0; i < strlen(duplex[1]); i++) {
                if (duplex[1].charAt(i) == 'A' || duplex[1].charAt(i) == 'T' ||
                        duplex[1].charAt(i) == 'C' || duplex[1].charAt(i) == 'G') {
                    replace(ret_str[1], i + 3, '|');
                } else {
                    replace(ret_str[1], i + 3, ' ');
                }
                truncate(ret_str[1], i + 4);
            }
            /* Clean Ends */
            ret_nr = strlen(ret_str[1]) - 1;
            while (ret_nr > 0 && ret_str[1].charAt(ret_nr) == ' ') {
                truncate(ret_str[1], ret_nr--);
            }
            /* Join bottom primer */
            strcpy(ret_str[2], "   ");
            strcat(ret_str[2], duplex[2]);
            ret_nr = 0;
            while (ret_nr < duplex[3].length()) {
                if (duplex[3].charAt(ret_nr) == 'A' || duplex[3].charAt(ret_nr) == 'T' ||
                        duplex[3].charAt(ret_nr) == 'C' || duplex[3].charAt(ret_nr) == 'G' ||
                        duplex[3].charAt(ret_nr) == '-') {
                    replace(ret_str[2], ret_nr + 3, duplex[3].charAt(ret_nr));
                }
                ret_nr++;
            }
            if (strlen(duplex[3]) > strlen(duplex[2])) {
                truncate(ret_str[2], strlen(duplex[3]) + 3);
            }
            /* Clean Ends */
            ret_nr = strlen(ret_str[2]) - 1;
            while (ret_nr > 0 && (ret_str[2].charAt(ret_nr) == ' ' || ret_str[2].charAt(ret_nr) == '-')) {
                truncate(ret_str[2], ret_nr--);
            }
            /* Write the 5' */
            ret_nr = 3;
            ret_pr_once = 1;
            while (ret_nr < ret_str[2].length() && ret_pr_once == 1) {
                if (ret_str[2].charAt(ret_nr) == 'A' || ret_str[2].charAt(ret_nr) == 'T' ||
                        ret_str[2].charAt(ret_nr) == 'C' || ret_str[2].charAt(ret_nr) == 'G' ||
                        ret_str[2].charAt(ret_nr) == '-') {
                    replace(ret_str[2], ret_nr - 3, '3');
                    replace(ret_str[2], ret_nr - 2, '\'');
                    ret_pr_once = 0;
                }
                ret_nr++;
            }

            append(ret_str[3], ret_para);
            append(ret_str[3], ret_str[0]);
            append(ret_str[3], " 3'\n");
            append(ret_str[3], ret_str[1]);
            append(ret_str[3], "\n");
            append(ret_str[3], ret_str[2]);
            append(ret_str[3], " 5'\n");

            return ret_str[3].toString();
        }

        return null;
    }

    /** prints ascii output of hairpin structure */
    String drawHairpin(int[] bp, double mh, double ms, final int mode, double temp, ThalResults o)
    {
        int ret_last_l, ret_first_r, ret_center, ret_left_end, ret_right_start, ret_left_len, ret_right_len;
        int ret_add_sp_l, ret_add_sp_r;
        char ret_center_char;
        String ret_para = "";
        /* Plain text */
        int i, N;
        N = 0;
        double mg, t;
        if (!isFinite(ms) || !isFinite(mh)) {
            o.temp = 0.0;
            o.msg = "Could not compute thermodynamic values";
            return "No predicted secondary structures for given sequences";
        } else {
            for (i = 1; i < len1; ++i) if (bp[i - 1] > 0) N++;
            // TODO: Shouldn't "(N / 2)" below be "(N / 2.0)" to avoid integer division errors?
            t = (mh / (ms + (((N / 2) - 1) * saltCorrection))) - ABSOLUTE_ZERO;
            mg = mh - (temp * (ms + (((N / 2) - 1) * saltCorrection)));
            ms = ms + (((N / 2) - 1) * saltCorrection);
            o.temp = t;
            o.delta_G = mg;
            o.delta_H = mh;
            o.delta_S = ms;
            if (mode == OUTPUT_MODE_FAST) return null;
            o.msg = String.format("Tm: %,12.2f \u00B0C\n\u0394G: %,12.2f cal/mol\n\u0394H: %,12.2f cal/mol\n\u0394S: %,12.2f cal/mol*K\n", t, mg, mh, ms);
        }
        /* plain-text output */
        char[] asciiRow = new char[len1];
        // for(i = 0; i < len1; ++i) asciiRow[i] = '0';
        for (i = 1; i < len1 + 1; ++i) {
            if (bp[i - 1] == 0) {
                asciiRow[(i - 1)] = '-';
            } else {
                if (bp[i - 1] > (i - 1)) {
                    asciiRow[(bp[i - 1] - 1)] = '\\';
                } else {
                    asciiRow[(bp[i - 1] - 1)] = '/';
                }
            }
        }
        if ((mode == OUTPUT_MODE_TEXT)) {
            StringBuilder sb = new StringBuilder("SEQ\t");
            for (i = 0; i < len1; ++i) sb.append(asciiRow[i]);
            sb.append(String.format("\nSTR\t%s\n", oligo1));
            return sb.toString();
        }
        if (mode == OUTPUT_MODE_GRAPHICAL) {
            StringBuilder ret_str = new StringBuilder();

            append(ret_str, ret_para);

            ret_last_l = -1;
            ret_first_r = -1;
            ret_center_char = '|';
            for (i = 0; i < len1; ++i) {
                if (asciiRow[i] == '/') {
                    ret_last_l = i;
                }
                if ((ret_first_r == -1) && (asciiRow[i] == '\\')) {
                    ret_first_r = i;
                }
            }
            ret_center = ret_first_r - ret_last_l;
            if (ret_center % 2 == 0) {
                /* ret_center is odd */
                ret_left_end = ret_last_l + (ret_first_r - ret_last_l) / 2 - 1;
                ret_center_char = (char) oligo1.charAt(ret_left_end + 1);
                ret_right_start = ret_left_end + 2;
            } else {
                /* ret_center is even */
                ret_left_end = ret_last_l + (ret_first_r - ret_last_l - 1) / 2;
                ret_right_start = ret_left_end + 1;
            }
            ret_left_len = ret_left_end + 1;
            ret_right_len = len1 - ret_right_start;
            ret_add_sp_l = 0;
            ret_add_sp_r = 0;
            if (ret_left_len > ret_right_len) {
                ret_add_sp_r = ret_left_len - ret_right_len + 1;
            }
            if (ret_right_len > ret_left_len) {
                ret_add_sp_l = ret_right_len - ret_left_len;
            }
            for (i = 0; i < ret_add_sp_l; i++) {
                append(ret_str, ' ');
            }
            append(ret_str, "5' ");
            for (i = 0; i < ret_left_len; i++) {
                append(ret_str, (char) oligo1.charAt(i));
            }
            append(ret_str, "\u2510\n   ");
            for (i = 0; i < ret_add_sp_l; i++) {
                append(ret_str, ' ');
            }
            for (i = 0; i < ret_left_len; i++) {
                if (asciiRow[i] == '/') {
                    append(ret_str, '|');
                } else {
                    append(ret_str, ' ');
                }
            }
            if (ret_center_char == '|') {
                append(ret_str, "\u2502");
            } else {
                append(ret_str, ret_center_char);
            }
            append(ret_str, "\n");
            for (i = 0; i < ret_add_sp_r - 1; i++) {
                append(ret_str, ' ');
            }
            append(ret_str, "3' ");
            for (i = len1 - 1; i > ret_right_start - 1; i--) {
                append(ret_str, (char) oligo1.charAt(i));
            }
            append(ret_str, "\u2518\n");

    /*
         append(ret_str, "SEQ ");
         for(i = 0; i < len1; ++i) {
           append(ret_str, asciiRow[i]);
         }
         append(ret_str, "\nSTR ");
         append(ret_str, (const char*) oligo1);
         append(ret_str, "\n");
    */
            return ret_str.toString();
        }
        return null;
    }

    private static void CHECK_ERROR(boolean COND, String MSG) {
        if (COND) throw new RuntimeException(MSG);
    }

    private static String LONG_SEQ_ERR_STR(int MAX_LEN) {
        return "Target sequence length > maximum allowed (" + MAX_LEN + ") in thermodynamic alignment";
    }

    /** Are the two numbers finite and very nearly equal (less than 1e-5 difference) ? */
    private static boolean equal(double a, double b) {
        if (!isFinite(a) || !isFinite(b))
            return false;
        return abs(a - b) < 1e-5;
    }

    /* 1 when numbers are nearly equal, 2 when unequal */
    // TODO: Should we use absolute value like equal() does?
    private static int DBL_EQ(double X, double Y) {
        return (X - Y < SMALL_NON_ZERO ? 1 : 2);
    }

    private static boolean isFinite(double x) {
        return Double.isFinite(x);
    }

    private static boolean isPositive(double x) {
        return x > 0;
    }

    private static void append(StringBuilder str, char c) {
        str.append(c);
    }

    private static void append(StringBuilder str, String c) {
        str.append(c);
    }

    private static void append(StringBuilder str, StringBuilder c) {
        str.append(c);
    }

    private static void replace(StringBuilder str, int i, char c) {
        str.replace(i, i + 1, Character.toString(c));
    }

    private static void truncate(StringBuilder str, int i) {
        str.delete(i, str.length());
    }

    private static void strcat(StringBuilder str, String c) {
        str.append(c);
    }

    private static void strcat(StringBuilder str, StringBuilder c) {
        str.append(c);
    }

    private static void strcpy(StringBuilder str, String c) {
        str.delete(0, str.length());
        str.append(c);
    }

    private static int strlen(StringBuilder str) {
        return str.length();
    }

    private static boolean empty_string(String a) {
        return a == null || a.length() == 0;
    }

    private static void printf(String fmt, Object... args) {
        System.out.printf(fmt, args);
    }

    private double EntropyDPT(int i, int j) {
        return entropyDPT[i][j];
    }

    private double EnthalpyDPT(int i, int j) {
        return enthalpyDPT[i][j];
    }

    private static void do_one(String left, String right, int mode) throws Exception {
        ThalArgs thal_arg = new ThalArgs();
        thal_arg.mode = mode;
        Thal thal = new Thal();
        if (right == null) {
            thal_arg.type = THAL_HAIRPIN;
            right = left;
        }
        else thal_arg.type = THAL_ANY;
        ThalResults r = thal.thal(left, right, thal_arg);

        if (!empty_string(r.msg)) System.out.printf("Message:\n%s", indent(r.msg));
        if (!empty_string(r.sec_struct)) System.out.printf("Secondary structure:\n%s", indent(r.sec_struct));
        if (mode == OUTPUT_MODE_FAST)
            System.out.printf("(Calculated values: tm %f, dG %f, dH %f, dS %f)\n",
                    r.temp, r.delta_G, r.delta_H, r.delta_S);
        System.out.print("\n");
    }

    private static String indent(String str) {
        String result = "";
        for (String line : str.split("\n")) {
            result += "    " + line + "\n";
        }
        return result;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 2) {
            System.out.print("Fast output mode:\n\n");
            do_one(args[0], args[1], OUTPUT_MODE_FAST);
            System.out.print("Simple output mode:\n\n");
            do_one(args[0], args[1], OUTPUT_MODE_TEXT);
            System.out.print("Graphical output mode:\n\n");
            do_one(args[0], args[1], OUTPUT_MODE_GRAPHICAL);
        } else if (args.length == 1) {
            System.out.print("Fast output mode:\n\n");
            do_one(args[0], null, OUTPUT_MODE_FAST);
            System.out.print("Simple output mode:\n\n");
            do_one(args[0], null, OUTPUT_MODE_TEXT);
            System.out.print("Graphical output mode:\n\n");
            do_one(args[0], null, OUTPUT_MODE_GRAPHICAL);
        } else {
            System.err.print("java -jar JavaPrimer3.jar oligo1 [oligo2]\n");
            System.exit(1);
        }
    }
}
