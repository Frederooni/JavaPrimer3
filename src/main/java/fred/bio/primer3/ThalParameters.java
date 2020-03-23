/*
 Copyright (c) 2018
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

/* This file is created by thal_parameters_c_create.pl
   Do not edit this file, edit the script instead!
 */

/* Converted to Java by Fred Long */

package fred.bio.primer3;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class ThalParameters {
    String dangle_dh;
    String dangle_ds;
    String loops_dh;
    String loops_ds;
    String stack_dh;
    String stack_ds;
    String stackmm_dh;
    String stackmm_ds;
    String tetraloop_dh;
    String tetraloop_ds;
    String triloop_dh;
    String triloop_ds;
    String tstack_tm_inf_ds;
    String tstack_dh;
    String tstack2_dh;
    String tstack2_ds;

    /**
     * Returns a ThalParameters instance that loads the thermodynamic values from given directory path.
     */
    public ThalParameters(String path) throws IOException {
        dangle_dh = readParamFile(path, "dangle.dh");
        dangle_ds = readParamFile(path, "dangle.ds");
        loops_dh = readParamFile(path, "loops.dh");
        loops_ds = readParamFile(path, "loops.ds");
        stack_dh = readParamFile(path, "stack.dh");
        stack_ds = readParamFile(path, "stack.ds");
        stackmm_dh = readParamFile(path, "stackmm.dh");
        stackmm_ds = readParamFile(path, "stackmm.ds");
        tetraloop_dh = readParamFile(path, "tetraloop.dh");
        tetraloop_ds = readParamFile(path, "tetraloop.ds");
        triloop_dh = readParamFile(path, "triloop.dh");
        triloop_ds = readParamFile(path, "triloop.ds");
        tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds");
        tstack_dh = readParamFile(path, "tstack.dh");
        tstack2_dh = readParamFile(path, "tstack2.dh");
        tstack2_ds = readParamFile(path, "tstack2.ds");
    }

    /** Returns a ThalParameters instance that loads the default thermodynamic values. */
    public ThalParameters() {
        try {
            Properties props = new Properties();
            InputStream is = ThalParameters.class.getResourceAsStream("/thal.properties");
            props.load(is);
            dangle_dh = props.getProperty("dangle_dh");
            dangle_ds = props.getProperty("dangle_ds");
            loops_dh = props.getProperty("loops_dh");
            loops_ds = props.getProperty("loops_ds");
            stack_dh = props.getProperty("stack_dh");
            stack_ds = props.getProperty("stack_ds");
            stackmm_dh = props.getProperty("stackmm_dh");
            stackmm_ds = props.getProperty("stackmm_ds");
            tetraloop_dh = props.getProperty("tetraloop_dh");
            tetraloop_ds = props.getProperty("tetraloop_ds");
            triloop_dh = props.getProperty("triloop_dh");
            triloop_ds = props.getProperty("triloop_ds");
            tstack_tm_inf_ds = props.getProperty("tstack_tm_inf_ds");
            tstack_dh = props.getProperty("tstack_dh");
            tstack2_dh = props.getProperty("tstack2_dh");
            tstack2_ds = props.getProperty("tstack2_ds");
        } catch (Exception e) {
            e.printStackTrace();
            System.err.printf(e + " - Could not load thal.properties, using built-in defaults.\n");
            set_default_thal_parameters(this);
        }
    }

    /**
     * Copy the default thermodynamic parameter strings to *a
     */
    public static void set_default_thal_parameters(ThalParameters a) {
        String dangle_dh = "0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n-500\n4700\n-4100\n-3800\n0\n" +
                "0\n0\n0\n0\n0\n0\n0\n-5900\n-2600\n-3200\n-5200\n0\n0\n0\n0\n0\n0\n" +
                "0\n0\n-2100\n-200\n-3900\n-4400\n0\n0\n0\n0\n0\n0\n0\n0\n-700\n4400\n-1600\n" +
                "2900\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n" +
                "0\n0\n0\n0\n0\n0\n0\n0\n-2900\n-4100\n-4200\n-200\n0\n0\n0\n0\n0\n" +
                "0\n0\n0\n-3700\n-4000\n-3900\n-4900\n0\n0\n0\n0\n0\n0\n0\n0\n-6300\n-4400\n" +
                "-5100\n-4000\n0\n0\n0\n0\n0\n0\n0\n0\n200\n600\n-1100\n-6900\n0\n0\n0\n" +
                "0\n0\n0\n0\n0\n0\n0\n0\n0\n";

        String dangle_ds = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-1.1\n14.2\n-13.1\n-12.6\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-16.5\n-7.4\n-10.4\n-15\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\n-3.9\n-0.1\n-11.2\n-13.1\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-0.8\n14.9\n-3.6\n" +
                "10.4\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7.6\n-13\n-15\n-0.5\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\n-10\n-11.9\n-10.9\n-13.8\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-17.1\n-12.6\n" +
                "-14\n-10.9\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n2.3\n3.3\n-1.6\n-20\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n";

        String loops_dh = "1\tinf\t0.0\tinf\n2\tinf\t0.0\tinf\n3\t0.0\t0.0\t0.0\n4\t0.0\t0.0\t0.0\n5\t0.0\t0.0\t0.0\n" +
                "6\t0.0\t0.0\t0.0\n7\t0.0\t0.0\t0.0\n8\t0.0\t0.0\t0.0\n9\t0.0\t0.0\t0.0\n10\t0.0\t0.0\t0.0\n" +
                "11\t0.0\t0.0\t0.0\n12\t0.0\t0.0\t0.0\n13\t0.0\t0.0\t0.0\n14\t0.0\t0.0\t0.0\n15\t0.0\t0.0\t0.0\n" +
                "16\t0.0\t0.0\t0.0\n17\t0.0\t0.0\t0.0\n18\t0.0\t0.0\t0.0\n19\t0.0\t0.0\t0.0\n20\t0.0\t0.0\t0.0\n" +
                "21\t0.0\t0.0\t0.0\n22\t0.0\t0.0\t0.0\n23\t0.0\t0.0\t0.0\n24\t0.0\t0.0\t0.0\n25\t0.0\t0.0\t0.0\n" +
                "26\t0.0\t0.0\t0.0\n27\t0.0\t0.0\t0.0\n28\t0.0\t0.0\t0.0\n29\t0.0\t0.0\t0.0\n30\t0.0\t0.0\t0.0\n" +
                "";

        String loops_ds = "1\t-1.0\t-12.89\t-1.0\n2\t-1.0\t-9.35\t-1.0\n3\t-10.31\t-9.99\t-11.28\n4\t-11.6\t-10.31\t-11.28\n5\t-12.89\t-10.64\t-10.64\n" +
                "6\t-14.18\t-11.28\t-12.89\n7\t-14.83\t-11.92\t-13.54\n8\t-15.47\t-12.57\t-13.86\n9\t-15.79\t-13.21\t-14.5\n10\t-15.79\t-13.86\t-14.83\n" +
                "11\t-16.26\t-14.32\t-15.29\n12\t-16.76\t-14.5\t-16.12\n13\t-17.15\t-14.89\t-16.5\n14\t-17.41\t-15.47\t-16.44\n15\t-17.74\t-15.81\t-16.77\n" +
                "16\t-18.05\t-16.12\t-17.08\n17\t-18.34\t-16.41\t-17.38\n18\t-18.7\t-16.76\t-17.73\n19\t-18.96\t-17.02\t-17.99\n20\t-19.02\t-17.08\t-18.37\n" +
                "21\t-19.25\t-17.32\t-18.61\n22\t-19.48\t-17.55\t-18.84\n23\t-19.7\t-17.76\t-19.05\n24\t-19.9\t-17.97\t-19.26\n25\t-20.31\t-18.05\t-19.66\n" +
                "26\t-20.5\t-18.24\t-19.85\n27\t-20.68\t-18.42\t-20.04\n28\t-20.86\t-18.6\t-20.21\n29\t-21.03\t-18.77\t-20.38\n30\t-21.28\t-19.02\t-20.31\n" +
                "";

        String stack_dh = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7900\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8400\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7800\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-8500\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\n-8000\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\n-10600\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\n-7800\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8200\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-9800\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8000\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-8400\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-7200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\n-8200\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\n-8500\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\n-7900\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\n";

        String stack_ds = "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.2\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.4\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-21.0\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-20.4\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.7\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\n-19.9\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\n-27.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\n-21.0\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.2\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-24.4\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-19.9\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-22.4\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n-21.3\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\ninf\ninf\n-22.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\ninf\ninf\n-22.7\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\ninf\n-22.2\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\ninf\n" +
                "inf\n";

        String stackmm_dh = "inf\ninf\ninf\n4700\ninf\ninf\ninf\n7600\ninf\ninf\ninf\n3000\n1200\n2300\n-600\ninf\ninf\n" +
                "inf\n-2900\ninf\ninf\ninf\n-700\ninf\ninf\ninf\n500\ninf\n5300\n-10\ninf\n700\ninf\n-900\n" +
                "inf\ninf\ninf\n600\ninf\ninf\ninf\n-4000\ninf\ninf\n-700\ninf\n-3100\n1000\n1200\ninf\ninf\n" +
                "inf\n5300\ninf\ninf\ninf\n-700\ninf\ninf\ninf\ninf\n-1200\n-2500\n-2700\ninf\ninf\ninf\n3400\n" +
                "inf\ninf\ninf\n6100\n-900\n1900\n-700\ninf\ninf\ninf\ninf\n1000\ninf\ninf\n5200\ninf\ninf\n" +
                "inf\n3600\ninf\n600\n-1500\ninf\n-800\ninf\ninf\n5200\ninf\ninf\n1900\ninf\ninf\ninf\n-1500\n" +
                "inf\ninf\n-4000\ninf\n-4900\n-4100\ninf\n-1500\ninf\ninf\n2300\ninf\ninf\ninf\n-10\ninf\ninf\n" +
                "inf\ninf\n-1500\n-2800\n-5000\n-1200\ninf\ninf\ninf\ninf\ninf\ninf\n700\n-2900\n5200\n-600\ninf\n" +
                "inf\ninf\ninf\n1600\ninf\ninf\ninf\n-1300\ninf\ninf\n-600\ninf\n-700\n3600\ninf\n2300\ninf\n" +
                "inf\n-6000\ninf\ninf\ninf\n-4400\ninf\ninf\n-700\ninf\ninf\n500\ninf\n-6000\n3300\ninf\n-4900\n" +
                "inf\ninf\ninf\n-2800\ninf\n5800\n-600\ninf\ninf\ninf\ninf\n5200\n-4400\n-2200\n-3100\ninf\ninf\n" +
                "inf\n-2500\ninf\n4100\ninf\ninf\n3400\n700\ninf\ninf\ninf\ninf\n1200\ninf\ninf\ninf\n-100\n" +
                "inf\ninf\ninf\n200\n7600\n6100\ninf\n1200\ninf\ninf\n2300\ninf\ninf\ninf\n3300\ninf\ninf\n" +
                "inf\n-2200\ninf\n3000\ninf\n1600\n-100\ninf\n-800\ninf\ninf\ninf\n-4100\ninf\n-1400\ninf\n-5000\n" +
                "inf\ninf\ninf\n1000\n-1300\n200\n700\ninf\ninf\ninf\n1000\ninf\n5800\ninf\n-2700\ninf\ninf\n" +
                "inf\n";

        String stackmm_ds = "inf\ninf\ninf\n12.9\ninf\ninf\ninf\n20.2\ninf\ninf\ninf\n7.4\n1.7\n4.6\n-2.3\ninf\ninf\n" +
                "inf\n-9.8\ninf\ninf\ninf\n-3.8\ninf\ninf\ninf\n3.2\ninf\n14.6\n-4.4\ninf\n0.2\ninf\n-4.2\n" +
                "inf\ninf\ninf\n-0.6\ninf\ninf\ninf\n-13.2\ninf\ninf\n-2.3\ninf\n-9.5\n0.9\n1.7\ninf\ninf\n" +
                "inf\n14.6\ninf\ninf\ninf\n-2.3\ninf\ninf\ninf\ninf\n-6.2\n-8.3\n-10.8\ninf\ninf\ninf\n8.0\n" +
                "inf\ninf\ninf\n16.4\n-4.2\n3.7\n-2.3\ninf\ninf\ninf\ninf\n0.7\ninf\ninf\n14.2\ninf\ninf\n" +
                "inf\n8.9\ninf\n-0.6\n-7.2\ninf\n-4.5\ninf\ninf\n13.5\ninf\ninf\n3.7\ninf\ninf\ninf\n-7.2\n" +
                "inf\ninf\n-13.2\ninf\n-15.3\n-11.7\ninf\n-6.1\ninf\ninf\n4.6\ninf\ninf\ninf\n-4.4\ninf\ninf\n" +
                "inf\ninf\n-6.1\n-8.0\n-15.8\n-6.2\ninf\ninf\ninf\ninf\ninf\ninf\n0.7\n-9.8\n14.2\n-1.0\ninf\n" +
                "inf\ninf\ninf\n3.6\ninf\ninf\ninf\n-5.3\ninf\ninf\n-1.0\ninf\n-3.8\n8.9\ninf\n5.4\ninf\n" +
                "inf\n-15.8\ninf\ninf\ninf\n-12.3\ninf\ninf\n-2.3\ninf\ninf\n3.2\ninf\n-15.8\n10.4\ninf\n-15.3\n" +
                "inf\ninf\ninf\n-8.0\ninf\n16.3\n-2.3\ninf\ninf\ninf\ninf\n13.5\n-12.3\n-8.4\n-9.5\ninf\ninf\n" +
                "inf\n-8.3\ninf\n9.5\ninf\n12.9\n8.0\n0.7\ninf\ninf\ninf\ninf\n0.7\ninf\ninf\ninf\n-1.7\n" +
                "inf\ninf\ninf\n-1.5\n20.2\n16.4\ninf\n0.7\ninf\ninf\n5.4\ninf\ninf\ninf\n10.4\ninf\ninf\n" +
                "inf\n-8.4\ninf\n7.4\ninf\n3.6\n-1.7\ninf\n-4.5\ninf\ninf\ninf\n-11.7\ninf\n-6.2\ninf\n-15.8\n" +
                "inf\ninf\ninf\n0.7\n-5.3\n-1.5\n0.2\ninf\ninf\ninf\n0.9\ninf\n16.3\ninf\n-10.8\ninf\ninf\n" +
                "inf\n";

        String tetraloop_dh = "AAAAAT\t500\nAAAACT\t700\nAAACAT\t1000\nACTTGT\t0\nAGAAAT\t-1100\n" +
                "AGAGAT\t-1100\nAGATAT\t-1500\nAGCAAT\t-1600\nAGCGAT\t-1100\nAGCTTT\t200\n" +
                "AGGAAT\t-1100\nAGGGAT\t-1100\nAGGGGT\t500\nAGTAAT\t-1600\nAGTGAT\t-1100\n" +
                "AGTTCT\t800\nATTCGT\t-200\nATTTGT\t0\nATTTTT\t-500\nCAAAAG\t500\n" +
                "CAAACG\t700\nCAACAG\t1000\nCAACCG\t0\nCCTTGG\t0\nCGAAAG\t-1100\n" +
                "CGAGAG\t-1100\nCGATAG\t-1500\nCGCAAG\t-1600\nCGCGAG\t-1100\nCGCTTG\t200\n" +
                "CGGAAG\t-1100\nCGGGAG\t-1000\nCGGGGG\t500\nCGTAAG\t-1600\nCGTGAG\t-1100\n" +
                "CGTTCG\t800\nCTTCGG\t-200\nCTTTGG\t0\nCTTTTG\t-500\nGAAAAC\t500\n" +
                "GAAACC\t700\nGAACAC\t1000\nGCTTGC\t0\nGGAAAC\t-1100\nGGAGAC\t-1100\n" +
                "GGATAC\t-1600\nGGCAAC\t-1600\nGGCGAC\t-1100\nGGCTTC\t200\nGGGAAC\t-1100\n" +
                "GGGGAC\t-1100\nGGGGGC\t500\nGGTAAC\t-1600\nGGTGAC\t-1100\nGGTTCC\t800\n" +
                "GTTCGC\t-200\nGTTTGC\t0\nGTTTTC\t-500\nTAAAAA\t500\nTAAACA\t700\n" +
                "TAACAA\t1000\nTCTTGA\t0\nTGAAAA\t-1100\nTGAGAA\t-1100\nTGATAA\t-1600\n" +
                "TGCAAA\t-1600\nTGCGAA\t-1100\nTGCTTA\t200\nTGGAAA\t-1100\nTGGGAA\t-1100\n" +
                "TGGGGA\t500\nTGTAAA\t-1600\nTGTGAA\t-1100\nTGTTCA\t800\nTTTCGA\t-200\n" +
                "TTTTGA\t0\nTTTTTA\t-500\n";

        String tetraloop_ds = "AAAAAT\t-650\nAAAACT\t1610\nAAACAT\t1610\nACTTGT\t4190\nAGAAAT\t1610\n" +
                "AGAGAT\t1610\nAGATAT\t1610\nAGCAAT\t1610\nAGCGAT\t1610\nAGCTTT\t1610\n" +
                "AGGAAT\t1610\nAGGGAT\t1610\nAGGGGT\t640\nAGTAAT\t1610\nAGTGAT\t1610\n" +
                "AGTTCT\t1610\nATTCGT\t1610\nATTTGT\t1610\nATTTTT\t1610\nCAAAAG\t-1290\n" +
                "CAAACG\t0\nCAACAG\t0\nCAACCG\t0\nCCTTGG\t2570\nCGAAAG\t0\n" +
                "CGAGAG\t0\nCGATAG\t0\nCGCAAG\t0\nCGCGAG\t0\nCGCTTG\t0\n" +
                "CGGAAG\t0\nCGGGAG\t0\nCGGGGG\t-970\nCGTAAG\t0\nCGTGAG\t0\n" +
                "CGTTCG\t0\nCTTCGG\t0\nCTTTGG\t0\nCTTTTG\t0\nGAAAAC\t-3230\n" +
                "GAAACC\t0\nGAACAC\t0\nGCTTGC\t2570\nGGAAAC\t0\nGGAGAC\t0\n" +
                "GGATAC\t0\nGGCAAC\t0\nGGCGAC\t0\nGGCTTC\t0\nGGGAAC\t0\n" +
                "GGGGAC\t0\nGGGGGC\t-970\nGGTAAC\t0\nGGTGAC\t0\nGGTTCC\t0\n" +
                "GTTCGC\t0\nGTTTGC\t0\nGTTTTC\t0\nTAAAAA\t320\nTAAACA\t1610\n" +
                "TAACAA\t1610\nTCTTGA\t4190\nTGAAAA\t1610\nTGAGAA\t1610\nTGATAA\t1610\n" +
                "TGCAAA\t1610\nTGCGAA\t1610\nTGCTTA\t1610\nTGGAAA\t1610\nTGGGAA\t1610\n" +
                "TGGGGA\t640\nTGTAAA\t1610\nTGTGAA\t1610\nTGTTCA\t1610\nTTTCGA\t1610\n" +
                "TTTTGA\t1610\nTTTTTA\t1610\n";

        String triloop_dh = "AGAAT\t-1500\nAGCAT\t-1500\nAGGAT\t-1500\nAGTAT\t-1500\nCGAAG\t-2000\n" +
                "CGCAG\t-2000\nCGGAG\t-2000\nCGTAG\t-2000\nGGAAC\t-2000\nGGCAC\t-2000\n" +
                "GGGAC\t-2000\nGGTAC\t-2000\nTGAAA\t-1500\nTGCAA\t-1500\nTGGAA\t-1500\n" +
                "TGTAA\t-1500\n";

        String triloop_ds = "AGAAT\t0\nAGCAT\t0\nAGGAT\t0\nAGTAT\t0\nCGAAG\t0\n" +
                "CGCAG\t0\nCGGAG\t0\nCGTAG\t0\nGGAAC\t0\nGGCAC\t0\n" +
                "GGGAC\t0\nGGTAC\t0\nTGAAA\t0\nTGCAA\t0\nTGGAA\t0\n" +
                "TGTAA\t0\n";

        String tstack_dh = "0\n0\n0\n-2500\n0\n0\n0\n-2700\n0\n0\n0\n-2400\n-3100\n-1600\n-1900\n0\n0\n" +
                "0\n-8000\n0\n0\n0\n-3200\n0\n0\n0\n-4600\n0\n-1800\n-100\n0\n-900\n0\n-4300\n" +
                "0\n0\n0\n-2700\n0\n0\n0\n-6000\n0\n0\n-2500\n0\n-1100\n-3200\n-3100\n0\n0\n" +
                "0\n-1800\n0\n0\n0\n-2500\n0\n0\n0\n0\n-2300\n-3500\n-2400\n0\n0\n0\n-2300\n" +
                "0\n0\n0\n-700\n-4300\n-2600\n-3900\n0\n0\n0\n0\n-700\n0\n0\n-5000\n0\n0\n" +
                "0\n-3900\n0\n-2700\n-2100\n0\n-3200\n0\n0\n-3000\n0\n0\n-2600\n0\n0\n0\n-2100\n" +
                "0\n0\n-6000\n0\n-3800\n-3800\n0\n-3900\n0\n0\n-1600\n0\n0\n0\n-100\n0\n0\n" +
                "0\n0\n-3900\n-6600\n-6100\n-2300\n0\n0\n0\n0\n0\n0\n-2000\n-8000\n-5000\n-4300\n0\n" +
                "0\n0\n0\n-1100\n0\n0\n0\n-3600\n0\n0\n-4300\n0\n-3200\n-3900\n0\n-4900\n0\n" +
                "0\n-700\n0\n0\n0\n-5900\n0\n0\n-3900\n0\n0\n-4600\n0\n-700\n-5700\n0\n-3800\n" +
                "0\n0\n0\n-6600\n0\n0\n-1900\n0\n0\n0\n0\n-3000\n-5900\n-7400\n-1100\n0\n0\n" +
                "0\n-3500\n0\n0\n0\n-2500\n-2300\n-2000\n-7200\n0\n0\n0\n-2500\n0\n0\n0\n-3900\n" +
                "0\n0\n0\n-3200\n-2700\n-700\n0\n-2500\n0\n0\n-4900\n0\n0\n0\n-5700\n0\n0\n" +
                "0\n-7400\n0\n-2400\n0\n-1100\n-3900\n0\n-3200\n0\n0\n0\n-3800\n0\n0\n0\n-6100\n" +
                "0\n0\n0\n-700\n-3600\n-3200\n-900\n0\n0\n0\n-3200\n0\n0\n0\n-2400\n0\n0\n" +
                "0\n";

        String tstack2_dh = "0\n0\n0\n-2500\n0\n0\n0\n-2700\n0\n0\n0\n-2400\n-3100\n-1600\n-1900\n-5000\n0\n" +
                "0\n-8000\n0\n0\n0\n-3200\n0\n0\n0\n-4600\n0\n-1800\n-100\n-6000\n-900\n0\n-4300\n" +
                "0\n0\n0\n-2700\n0\n0\n0\n-6000\n0\n0\n-2500\n-6000\n-1100\n-3200\n-3100\n0\n0\n" +
                "0\n-1800\n0\n0\n0\n-2500\n0\n0\n0\n-5000\n-2300\n-3500\n-2400\n0\n0\n0\n-2300\n" +
                "0\n0\n0\n-700\n-4300\n-2600\n-3900\n-6000\n0\n0\n0\n-700\n0\n0\n-5000\n0\n0\n" +
                "0\n-3900\n0\n-2700\n-2100\n-7000\n-3200\n0\n0\n-3000\n0\n0\n-2600\n0\n0\n0\n-2100\n" +
                "0\n0\n-6000\n-7000\n-3800\n-3800\n0\n-3900\n0\n0\n-1600\n0\n0\n0\n-100\n0\n0\n" +
                "0\n-6000\n-3900\n-6600\n-6100\n-2300\n0\n0\n0\n0\n0\n0\n-2000\n-8000\n-5000\n-4300\n-6000\n" +
                "0\n0\n0\n-1100\n0\n0\n0\n-3600\n0\n0\n-4300\n0\n-3200\n-3900\n-7000\n-4900\n0\n" +
                "0\n-700\n0\n0\n0\n-5900\n0\n0\n-3900\n0\n0\n-4600\n-7000\n-700\n-5700\n0\n-3800\n" +
                "0\n0\n0\n-6600\n0\n0\n-1900\n0\n0\n0\n-6000\n-3000\n-5900\n-7400\n-1100\n0\n0\n" +
                "0\n-3500\n0\n0\n0\n-2500\n-2300\n-2000\n-5000\n0\n0\n0\n-2500\n0\n0\n0\n-3900\n" +
                "0\n0\n0\n-3200\n-2700\n-700\n-6000\n-2500\n0\n0\n-4900\n0\n0\n0\n-5700\n0\n0\n" +
                "0\n-7400\n0\n-2400\n-6000\n-1100\n-3900\n0\n-3200\n0\n0\n0\n-3800\n0\n0\n0\n-6100\n" +
                "0\n0\n-5000\n-700\n-3600\n-3200\n-900\n0\n0\n0\n-3200\n0\n0\n0\n-2400\n0\n0\n" +
                "0\n";

        String tstack2_ds = "inf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-7.0\ninf\ninf\ninf\n-5.8\n-7.8\n-4.0\n-4.4\n-13.5\ninf\n" +
                "inf\n-22.5\ninf\ninf\ninf\n-7.1\ninf\ninf\ninf\n-11.4\ninf\n-3.8\n-0.5\n-16.1\n-1.7\ninf\n-10.7\n" +
                "inf\ninf\ninf\n-6.0\ninf\ninf\ninf\n-15.5\ninf\ninf\n-5.9\n-16.1\n-2.1\n-8.7\n-7.8\ninf\ninf\n" +
                "inf\n-3.8\ninf\ninf\ninf\n-5.9\ninf\ninf\ninf\n-13.6\n-6.3\n-9.4\n-6.5\ninf\ninf\ninf\n-5.9\n" +
                "inf\ninf\ninf\n-1.3\n-10.7\n-5.9\n-9.6\n-16.1\ninf\ninf\ninf\n-1.2\ninf\ninf\n-13.8\ninf\ninf\n" +
                "inf\n-10.6\ninf\n-6.0\n-5.1\n-19.3\n-8.0\ninf\ninf\n-7.8\ninf\ninf\n-5.9\ninf\ninf\ninf\n-5.1\n" +
                "inf\ninf\n-15.5\n-19.3\n-9.5\n-9.0\ninf\n-10.6\ninf\ninf\n-4.0\ninf\ninf\ninf\n-0.5\ninf\ninf\n" +
                "inf\n-16.1\n-10.6\n-18.7\n-16.9\n-6.3\ninf\ninf\ninf\ninf\ninf\ninf\n-4.7\n-22.5\n-13.8\n-11.1\n-16.1\n" +
                "inf\ninf\ninf\n-2.7\ninf\ninf\ninf\n-9.8\ninf\ninf\n-11.1\ninf\n-7.1\n-10.6\n-19.3\n-13.5\ninf\n" +
                "inf\n-19.2\ninf\ninf\ninf\n-16.1\ninf\ninf\n-9.6\ninf\ninf\n-11.4\n-19.3\n-19.2\n-15.9\ninf\n-9.5\n" +
                "inf\ninf\ninf\n-18.7\ninf\ninf\n-4.4\ninf\ninf\ninf\n-16.1\n-7.8\n-16.1\n-21.2\n-2.1\ninf\ninf\n" +
                "inf\n-9.4\ninf\ninf\ninf\n-6.3\n-5.9\n-4.7\n-14.2\ninf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-10.5\n" +
                "inf\ninf\ninf\n-8.9\n-7.0\n-1.3\n-16.1\n-6.3\ninf\ninf\n-13.5\ninf\ninf\ninf\n-15.9\ninf\ninf\n" +
                "inf\n-21.2\ninf\n-5.8\n-16.1\n-2.7\n-10.5\ninf\n-8.0\ninf\ninf\ninf\n-9.0\ninf\ninf\ninf\n-16.9\n" +
                "inf\ninf\n-13.5\n-1.2\n-9.8\n-8.9\n-1.7\ninf\ninf\ninf\n-8.7\ninf\ninf\ninf\n-6.5\ninf\ninf\n" +
                "inf\n";

        String tstack_tm_inf_ds = "inf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-7.0\ninf\ninf\ninf\n-5.8\n-7.8\n-4.0\n-4.4\ninf\ninf\n" +
                "inf\n-22.5\ninf\ninf\ninf\n-7.1\ninf\ninf\ninf\n-11.4\ninf\n-3.8\n-0.5\ninf\n-1.7\ninf\n-10.7\n" +
                "inf\ninf\ninf\n-6.0\ninf\ninf\ninf\n-15.5\ninf\ninf\n-5.9\ninf\n-2.1\n-8.7\n-7.8\ninf\ninf\n" +
                "inf\n-3.8\ninf\ninf\ninf\n-5.9\ninf\ninf\ninf\ninf\n-6.3\n-9.4\n-6.5\ninf\ninf\ninf\n-5.9\n" +
                "inf\ninf\ninf\n-1.3\n-10.7\n-5.9\n-9.6\ninf\ninf\ninf\ninf\n-1.2\ninf\ninf\n-13.8\ninf\ninf\n" +
                "inf\n-10.6\ninf\n-6.0\n-5.1\ninf\n-8.0\ninf\ninf\n-7.8\ninf\ninf\n-5.9\ninf\ninf\ninf\n-5.1\n" +
                "inf\ninf\n-15.5\ninf\n-9.5\n-9.0\ninf\n-10.6\ninf\ninf\n-4.0\ninf\ninf\ninf\n-0.5\ninf\ninf\n" +
                "inf\ninf\n-10.6\n-18.7\n-16.9\n-6.3\ninf\ninf\ninf\ninf\ninf\ninf\n-4.7\n-22.5\n-13.8\n-11.1\ninf\n" +
                "inf\ninf\ninf\n-2.7\ninf\ninf\ninf\n-9.8\ninf\ninf\n-11.1\ninf\n-7.1\n-10.6\ninf\n-13.5\ninf\n" +
                "inf\n-19.2\ninf\ninf\ninf\n-16.1\ninf\ninf\n-9.6\ninf\ninf\n-11.4\ninf\n-19.2\n-15.9\ninf\n-9.5\n" +
                "inf\ninf\ninf\n-18.7\ninf\ninf\n-4.4\ninf\ninf\ninf\ninf\n-7.8\n-16.1\n-21.2\n-2.1\ninf\ninf\n" +
                "inf\n-9.4\ninf\ninf\ninf\n-6.3\n-5.9\n-4.7\ninf\ninf\ninf\ninf\n-6.3\ninf\ninf\ninf\n-10.5\n" +
                "inf\ninf\ninf\n-8.9\n-7.0\n-1.3\ninf\n-6.3\ninf\ninf\n-13.5\ninf\ninf\ninf\n-15.9\ninf\ninf\n" +
                "inf\n-21.2\ninf\n-5.8\ninf\n-2.7\n-10.5\ninf\n-8.0\ninf\ninf\ninf\n-9.0\ninf\ninf\ninf\n-16.9\n" +
                "inf\ninf\ninf\n-1.2\n-9.8\n-8.9\n-1.7\ninf\ninf\ninf\n-8.7\ninf\ninf\ninf\n-6.5\ninf\ninf\n" +
                "inf\n";

        a.dangle_dh = dangle_dh;
        a.dangle_ds = dangle_ds;
        a.loops_dh = loops_dh;
        a.loops_ds = loops_ds;
        a.stack_dh = stack_dh;
        a.stack_ds = stack_ds;
        a.stackmm_dh = stackmm_dh;
        a.stackmm_ds = stackmm_ds;
        a.tetraloop_dh = tetraloop_dh;
        a.tetraloop_ds = tetraloop_ds;
        a.triloop_dh = triloop_dh;
        a.triloop_ds = triloop_ds;
        a.tstack_tm_inf_ds = tstack_tm_inf_ds;
        a.tstack_dh = tstack_dh;
        a.tstack2_dh = tstack2_dh;
        a.tstack2_ds = tstack2_ds;
    }

    /** Return the contents of the parameter file as a string. */
    private static String readParamFile(final String dirname, final String fname) throws IOException {
        FileReader is = new FileReader(new File(dirname, fname));
        char[] buffer = new char[256];
        int c;
        StringBuilder sb = new StringBuilder();
        while ((c = is.read(buffer)) >= 0) sb.append(buffer, 0, c);
        return sb.toString();
    }

    private static void exportProperties(String name, String value) {
        System.out.printf("%s=%s\n", name, value.replaceAll("\\n", "\\\\n"));
    }

    public static void main(String[] args) {
        ThalParameters p = new ThalParameters();
        exportProperties("dangle_dh", p.dangle_dh);
        exportProperties("dangle_ds", p.dangle_ds);
        exportProperties("loops_dh", p.loops_dh);
        exportProperties("loops_ds", p.loops_ds);
        exportProperties("stack_dh", p.stack_dh);
        exportProperties("stack_ds", p.stack_ds);
        exportProperties("stackmm_dh", p.stackmm_dh);
        exportProperties("stackmm_ds", p.stackmm_ds);
        exportProperties("tetraloop_dh", p.tetraloop_dh);
        exportProperties("tetraloop_ds", p.tetraloop_ds);
        exportProperties("triloop_dh", p.triloop_dh);
        exportProperties("triloop_ds", p.triloop_ds);
        exportProperties("tstack_tm_inf_ds", p.tstack_tm_inf_ds);
        exportProperties("tstack_dh", p.tstack_dh);
        exportProperties("tstack2_dh", p.tstack2_dh);
        exportProperties("tstack2_ds", p.tstack2_ds);
    }
}
