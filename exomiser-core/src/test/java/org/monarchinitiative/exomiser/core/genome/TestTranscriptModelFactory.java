/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2019 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.core.genome;

import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.hgnc.AltGeneIDType;
import de.charite.compbio.jannovar.reference.HG19RefDictBuilder;
import de.charite.compbio.jannovar.reference.Strand;
import de.charite.compbio.jannovar.reference.TranscriptModel;
import org.monarchinitiative.exomiser.core.model.GeneIdentifier;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

import static de.charite.compbio.jannovar.reference.Strand.FWD;
import static de.charite.compbio.jannovar.reference.Strand.REV;

/**
 * Allows the easy creation of transcript models from the UCSC knownGenes.txt.gz lines.
 * The sequences are the output of the 5'UTR Exons, CDS Exons and 3'UTR Exons from UCSC
 *
 * @author Manuel Holtgrewe <manuel.holtgrewe@charite.de>
 * @author Jules Jacobsen <julesjacobsen@sanger.ac.uk>
 */
class TestTranscriptModelFactory {

    private static final Logger logger = LoggerFactory.getLogger(TestTranscriptModelFactory.class);

    private static final ReferenceDictionary REF_DICT = HG19RefDictBuilder.build();

    /**
     * @return {@link TranscriptModel} for gene FGFR2.
     */
    public static TranscriptModel buildTMForFGFR2() {
        return new TranscriptModelBuilder()
                .geneIdentifier(TestGeneFactory.FGFR2_IDENTIFIER)
                .knownGeneLine("uc021pzz.1\tchr10\t-\t123237843\t123357972\t123239370\t123353331\t18\t123237843,123243211,123244908,123246867,123247504,123256045,123258008,123260339,123263303,123274630,123276832,123279492,123298105,123310803,123324015,123324951,123353222,123357475,\t123239535,123243317,123245046,123246938,123247627,123256236,123258119,123260461,123263455,123274833,123276977,123279683,123298229,123310973,123324093,123325218,123353481,123357972,\tP21802\tuc021pzz.1")
                .mRnaSequence("ggcggcggctggaggagagcgcggtggagagccgagcgggcgggcggcgggtgcggagcgggcgagggagcgcgcgcggccgccacaaagctcgggcgccgcggggctgcatgcggcgtacctggcccggcgcggcgactgctctccgggctggcgggggccggccgcgagccccgggggccccgaggccgcagcttgcctgcgcgctctgagccttcgcaactcgcgagcaaagtttggtggaggcaacgccaagcctgagtcctttcttcctctcgttccccaaatccgagggcagcccgcgggcgtcatgcccgcgctcctccgcagcctggggtacgcgtgaagcccgggaggcttggcgccggcgaagacccaaggaccactcttctgcgtttggagttgctccccgcaaccccgggctcgtcgctttctccatcccgacccacgcggggcgcggggacaacacaggtcgcggaggagcgttgccattcaagtgactgcagcagcagcggcagcgcctcggttcctgagcccaccgcaggctgaaggcattgcgcgtagtccatgcccgtagaggaagtgtgcagatgggattaacgtccacatggagatatggaagaggaccggggattggtaccgtaaccatggtcagctggggtcgtttcatctgcctggtcgtggtcaccatggcaaccttgtccctggcccggccctccttcagtttagttgaggataccacattagagccagaagagccaccaaccaaataccaaatctctcaaccagaagtgtacgtggctgcgccaggggagtcgctagaggtgcgctgcctgttgaaagatgccgccgtgatcagttggactaaggatggggtgcacttggggcccaacaataggacagtgcttattggggagtacttgcagataaagggcgccacgcctagagactccggcctctatgcttgtactgccagtaggactgtagacagtgaaacttggtacttcatggtgaatgtcacagatgccatctcatccggagatgatgaggatgacaccgatggtgcggaagattttgtcagtgagaacagtaacaacaagagagcaccatactggaccaacacagaaaagatggaaaagcggctccatgctgtgcctgcggccaacactgtcaagtttcgctgcccagccggggggaacccaatgccaaccatgcggtggctgaaaaacgggaaggagtttaagcaggagcatcgcattggaggctacaaggtacgaaaccagcactggagcctcattatggaaagtgtggtcccatctgacaagggaaattatacctgtgtagtggagaatgaatacgggtccatcaatcacacgtaccacctggatgttgtggagcgatcgcctcaccggcccatcctccaagccggactgccggcaaatgcctccacagtggtcggaggagacgtagagtttgtctgcaaggtttacagtgatgcccagccccacatccagtggatcaagcacgtggaaaagaacggcagtaaatacgggcccgacgggctgccctacctcaaggttctcaaggccgccggtgttaacaccacggacaaagagattgaggttctctatattcggaatgtaacttttgaggacgctggggaatatacgtgcttggcgggtaattctattgggatatcctttcactctgcatggttgacagttctgccagcgcctggaagagaaaaggagattacagcttccccagactacctggagatagccatttactgcataggggtcttcttaatcgcctgtatggtggtaacagtcatcctgtgccgaatgaagaacacgaccaagaagccagacttcagcagccagccggctgtgcacaagctgaccaaacgtatccccctgcggagacaggtaacagtttcggctgagtccagctcctccatgaactccaacaccccgctggtgaggataacaacacgcctctcttcaacggcagacacccccatgctggcaggggtctccgagtatgaacttccagaggacccaaaatgggagtttccaagagataagctgacactgggcaagcccctgggagaaggttgctttgggcaagtggtcatggcggaagcagtgggaattgacaaagacaagcccaaggaggcggtcaccgtggccgtgaagatgttgaaagatgatgccacagagaaagacctttctgatctggtgtcagagatggagatgatgaagatgattgggaaacacaagaatatcataaatcttcttggagcctgcacacaggatgggcctctctatgtcatagttgagtatgcctctaaaggcaacctccgagaatacctccgagcccggaggccacccgggatggagtactcctatgacattaaccgtgttcctgaggagcagatgaccttcaaggacttggtgtcatgcacctaccagctggccagaggcatggagtacttggcttcccaaaaatgtattcatcgagatttagcagccagaaatgttttggtaacagaaaacaatgtgatgaaaatagcagactttggactcgccagagatatcaacaatatagactattacaaaaagaccaccaatgggcggcttccagtcaagtggatggctccagaagccctgtttgatagagtatacactcatcagagtgatgtctggtccttcggggtgttaatgtgggagatcttcactttagggggctcgccctacccagggattcccgtggaggaactttttaagctgctgaaggaaggacacagaatggataagccagccaactgcaccaacgaactgtacatgatgatgagggactgttggcatgcagtgccctcccagagaccaacgttcaagcagttggtagaagacttggatcgaattctcactctcacaaccaatgaggaatacttggacctcagccaacctctcgaacagtattcacctagttaccctgacacaagaagttcttgttcttcaggagatgattctgttttttctccagaccccatgccttacgaaccatgccttcctcagtatccacacataaacggcagtgttaaaacatgaatgactgtgtctgcctgtccccaaacaggacagcactgggaacctagctacactgagcagggagaccatgcctcccagagcttgttgtctccacttgtatatatggatcagaggagtaaataattggaaaagtaatcagcatatgtgtaaagatttatacagttgaaaacttgtaatcttccccaggaggagaagaaggtttctggagcagtggactgccacaagccaccatgtaacccctctcacctgccgtgcgtactggctgtggaccagtaggactcaaggtggacgtgcgttctgccttccttgttaattttgtaataattggagaagatttatgtcagcacacacttacagagcacaaatgcagtatataggtgctggatgtatgtaaatatattcaaattatgtataaatatatattatatatttacaaggagttattttttgtattgattttaaatggatgtcccaatgcacctagaaaattggtctctctttttttaatagctatttgctaaatgctgttcttacacataatttcttaattttcaccgagcagaggtggaaaaatacttttgctttcagggaaaatggtataacgttaatttattaataaattggtaatatacaaaacaattaatcatttatagttttttttgtaatttaagtggcatttctatgcaggcagcacagcagactagttaatctattgcttggacttaactagttatcagatcctttgaaaagagaatatttacaatatatgactaatttggggaaaatgaagttttgatttatttgtgtttaaatgctgctgtcagacgattgttcttagacctcctaaatgccccatattaaaagaactcattcataggaaggtgtttcattttggtgtgcaaccctgtcattacgtcaacgcaacgtctaactggacttcccaagataaatggtaccagcgtcctcttaaaagatgccttaatccattccttgaggacagaccttagttgaaatgatagcagaatgtgcttctctctggcagctggccttctgcttctgagttgcacattaatcagattagcctgtattctcttcagtgaattttgataatggcttccagactctttggcgttggagacgcctgttaggatcttcaagtcccatcatagaaaattgaaacacagagttgttctgctgatagttttggggatacgtccatctttttaagggattgctttcatctaattctggcaggacctcaccaaaagatccagcctcatacctacatcagacaaaatatcgccgttgttccttctgtactaaagtattgtgttttgctttggaaacacccactcactttgcaatagccgtgcaagatgaatgcagattacactgatcttatgtgttacaaaattggagaaagtatttaataaaacctgttaatttttatactgacaataaaaatgtttctacagatattaatgttaacaagacaaaataaatgtcacgcaacttatttttttaataaaaaaaaaaaaaaa")
                .build();
    }

    /**
     * @return {@link TranscriptModel} for gene SHH.
     */
    public static TranscriptModel buildTMForSHH() {
        return new TranscriptModelBuilder()
                .geneIdentifier(TestGeneFactory.SHH_IDENTIFIER)
                .knownGeneLine("uc003wmk.1\tchr7\t-\t155595557\t155604967\t155595593\t155604816\t3\t155595557,155598989,155604516,\t155596420,155599251,155604967,\tQ15465\tuc003wmk.1")
                .mRnaSequence("gcgaggcagccagcgagggagagagcgagcgggcgagccggagcgaggaagggaaagcgcaagagagagcgcacacgcacacacccgccgcgcgcactcgcgcacggacccgcacggggacagctcggaagtcatcagttccatgggcgagatgctgctgctggcgagatgtctgctgctagtcctcgtctcctcgctgctggtatgctcgggactggcgtgcggaccgggcagggggttcgggaagaggaggcaccccaaaaagctgacccctttagcctacaagcagtttatccccaatgtggccgagaagaccctaggcgccagcggaaggtatgaagggaagatctccagaaactccgagcgatttaaggaactcacccccaattacaaccccgacatcatatttaaggatgaagaaaacaccggagcggacaggctgatgactcagaggtgtaaggacaagttgaacgctttggccatctcggtgatgaaccagtggccaggagtgaaactgcgggtgaccgagggctgggacgaagatggccaccactcagaggagtctctgcactacgagggccgcgcagtggacatcaccacgtctgaccgcgaccgcagcaagtacggcatgctggcccgcctggcggtggaggccggcttcgactgggtgtactacgagtccaaggcacatatccactgctcggtgaaagcagagaactcggtggcggccaaatcgggaggctgcttcccgggctcggccacggtgcacctggagcagggcggcaccaagctggtgaaggacctgagccccggggaccgcgtgctggcggcggacgaccagggccggctgctctacagcgacttcctcactttcctggaccgcgacgacggcgccaagaaggtcttctacgtgatcgagacgcgggagccgcgcgagcgcctgctgctcaccgccgcgcacctgctctttgtggcgccgcacaacgactcggccaccggggagcccgaggcgtcctcgggctcggggccgccttccgggggcgcactggggcctcgggcgctgttcgccagccgcgtgcgcccgggccagcgcgtgtacgtggtggccgagcgtgacggggaccgccggctcctgcccgccgctgtgcacagcgtgaccctaagcgaggaggccgcgggcgcctacgcgccgctcacggcccagggcaccattctcatcaaccgggtgctggcctcgtgctacgcggtcatcgaggagcacagctgggcgcaccgggccttcgcgcccttccgcctggcgcacgcgctcctggctgcactggcgcccgcgcgcacggaccgcggcggggacagcggcggcggggaccgcgggggcggcggcggcagagtagccctaaccgctccaggtgctgccgacgctccgggtgcgggggccaccgcgggcatccactggtactcgcagctgctctaccaaataggcacctggctcctggacagcgaggccctgcacccgctgggcatggcggtcaagtccagctgaagccggggggccgggggaggggcgcgggagggggcg")
                .build();
    }

    /**
     * GNRHR2 overlaps with RBM8A.
     *
     * http://genome-euro.ucsc.edu/cgi-bin/hgc?hgsid=210987005_tfjW2xHSIt6ihH2y9UG6IxbvBx4B&g=htcGeneInGenome&i=uc009wiv.3&c=chr1&l=145509751&r=145515899&o=knownGene&table=knownGene
     *
     * @return {@link TranscriptModel} for gene GNRHR2.
     */
    public static TranscriptModel buildTMForGNRHR2A() {
        return new TranscriptModelBuilder()
                .geneIdentifier(TestGeneFactory.GNRHR2_IDENTIFIER)
                .knownGeneLine("uc009wiv.3\tchr1\t-\t145509751\t145515899\t145509859\t145510794\t4\t145509751,145510727,145515188,145515669,\t145510278,145510938,145515409,145515899,\tQ96P88\tuc009wiv.3")
                .mRnaSequence("atcttcgttccctgcagaaccttgacagttgaacaagtgacctcctccagaacagatggagagtctccagaagcagaggctttagtgaacgaaattcgcaataatcagctccagatcctgaaaaggagggcgaagaatcagtggccaaagctaaccgcttcatacccacacttcatcctcctcagtttctctccaggccaccatgtctgcaggcaacggcaccccttgggatgccacctggaatatcactgttcaatggctggctgtggacatcgcatgtcggacactgatgttcctgaaactaatggccacgtattctgcagctttcctgcctgtggtcattggattggaccgccaggcagcagtactcaacccgcttggatcccgttcaggtgtaaggaaacttctgggggcagcctggggacttagtttcctgcttgccttcccccagctgttcctgttccacacggtccactgagctggcccagtccctttcactcagtgtgtcaccaaaggcagcttcaaggctcaatggcaagagaccacctataacctcttcaccttctgctgcctctttctgctgccactgactgccatggccatctgctatagccgcattgtcctcagtgtgtccaggccccagacaaggaaggggagccatgcccctgctggtgaatttgccctcccccgctcctttgacaattgtccccgtgttcgtctccgggccctgagactggccctgcttatcttgctgaccttcatcctctgctggacaccttattacctactgggtatgtggtactggttctcccccaccatgctaactgaagtccctcccagcctgagccacatccttttcctcttgggcctcctcaatgctcctttggatcctctcctctatggggccttcacccttggctgccgaagagggcaccaagaacttagtatagactcttctaaagaagggtctgggagaatgctccaagaggagattcatgcctttagacagctggaagtacaaaaaactgtgacatcaagaagggcaggagaaacaaaaggcatttctataacatctatctgatcctaacagagtatgtaggaacagaatagtaagtctttagtgccataagatcttaacatctcacttctactcctgctctcctagttccccccaaaaaagaaatactga")
                .build();
    }

    /**
     * RBM8A overlaps with GNRHR2.
     *
     * @return {@link TranscriptModel} for gene RBM8A.
     */
    public static TranscriptModel buildTMForRBM8A() {
        return new TranscriptModelBuilder()
                .geneIdentifier(TestGeneFactory.RBM8A_IDENTIFIER)
                .knownGeneLine("uc001ent.2\tchr1\t+\t145507556\t145513535\t145507666\t145509211\t6\t145507556,145508015,145508206,145508474,145508915,145509165,\t145507733,145508075,145508284,145508611,145509052,145513535,\tQ9Y5S9\tuc001ent.2")
                .mRnaSequence("agagttagcctttgattggtcagcttgactggcgacctttcccctctgcgacagtttcccgaggtacctagtgtctgagcggcacagacgagatctcgatcgaaggcgagatggcggacgtgctagatcttcacgaggctgggggcgaagatttcgccatggatgaggatggggacgagagcattcacaaactgaaagaaaaagcgaagaaacggaagggtcgcggctttggctccgaagaggggtcccgagcgcggatgcgtgaggattatgacagcgtggagcaggatggcgatgaacccggaccacaacgctctgttgaaggctggattctctttgtaactggagtccatgaggaagccaccgaagaagacatacacgacaaattcgcagaatatggggaaattaaaaacattcatctcaacctcgacaggcgaacaggatatctgaaggggtatactctagttgaatatgaaacatacaaggaagcccaggctgctatggagggactcaatggccaggatttgatgggacagcccatcagcgttgactggtgttttgttcggggtccaccaaaaggcaagaggagaggtggccgaagacgcagcagaagtccagaccggagacgtcgctgacaggtcctctgttgtccaggtgttctcttcaagattccatttgaccatgcagccttggacaaataggactggggtggaacttgctgtgtttatatttaatctcttaccgtatatgcgtagtatttgagttgcgaataaatgttccatttttgttttctacatttaatgttactttcctgtcctaaaattgaaagttctaaagcatagcaaggctgtatggatcattgtgaagatacttctagggactgaactctatgtatttcttttttttcttttttttgagatagagtcttgctgtgttacccagggtggattgcagctgatcatagctcactgcagcttcaaactcttgggctcaagccatccttctgcctcactgtccctagtagttgggattacaggcacatgccaccatgcccagctaaatttttaatatttttgtagagatggggtcttgctgtgttacctgggctagttatgtgagtttctatattagacatagtctcaagtttcaggtagggtttaaagtagagacactggtcagtatttcttttttggggggaactaggagagcaggagtagaagtgagatgttaagatcttatggcactaaagacttactattctgttcctacatactctgttaggatcagatagatgttatagaaatgccttttgtttctcctgcccttcttgatgtcacagttttttgtacttccagctgtctaaaggcatgaatctcctcttggagcattctcccagacccttctttagaagagtctatactaagttcttggtgccctcttcggcagccaagggtgaaggccccatagaggagaggatccaaaggagcattgaggaggcccaagaggaaaaggatgtggctcaggctgggagggacttcagttagcatggtgggggagaaccagtaccacatacccagtaggtaataaggtgtccagcagaggatgaaggtcagcaagataagcagggccagtctcagggcccggagacgaacacggggacaattgtcaaaggagcgggggagggcaaattcaccagcaggggctaggaatttagaaaatatactgtaattcagacactcagcttctgatctgagtatagggtgaattgatggaggggcatagctagtgagacagagctcacctcctacaaggaggagaatgttgcaaaccgttttccccttcccaacctgggactatatgatttcttacccccagggattatgatagaaatatgaagccaccaagtctagacttgatggtgttcaagaataaataatactgattgcctccctagtccttgtccagctaactcagctgtttataattgaagggattcaacaaaattatctctagcatcaggtgctagacatggttagaatctcaccatggtttagtgactggtagatagctattaggtaggtagataaataaatgatgctagaggcaacaggtctagggttaaggattaaggcctgggaattggagtctcaccatggctccccttccttgtctggggcctggacacactgaggacaatgcggctatagcagatggccatggcagtcagtggcagcagaaagaggcagcagaaggtgaagaggttataggtggtctcttgccattgagccttgaagctgcctttggtgacacactgagtgaaagggactgggccagctcagtggaccgtgtggaacaggaacagctggagtggagttgaggactattagaactggttcccctcaccacccaacctacccacctatgtcatactgtctcctcccaattcatccttaattccaagtgaagcagcacagtgctgagaaacagttcatccatggtgccatgttaaagaagttggaaatatatcttgaaaatcctatcttccttttaggcttgaatatgatgctgaacagtaagtttgttaaatcttggaacttaaaacaatcctgctttctcaagtactattctaacattgcgctttataagggatgatatttctaccacctcactcatatttttagctgaaatgattttcctggtatgtctgttattttgtggaaaaagaaatattgtgtaaaatgggtgctgccaaaattccaggccattttgcagggactctgaagtgacctttagtagtaatagtcttatgtgcagtaactataatggtaaagaatgttaaataataaaatttaacattttccaaatgctattgggctgcccctccccctttttgttaaattgctgggttttccaactgaatcagtaaaaactatttctgtttagagctacaaggttaaagtgcctgctttccagtaatggagattgagtcactattaatttgataaaaggtaagctcagtaggcatcagattcctagatacaaggcatttgggaaagtgattttagcagacatgagggacatttaggaaagatgaatagtttcagcctaagagaattttgtgaactgtttggagttacgatcaggctactctgagctagttgggaaatggtctttcctcttcccatctcttgcattcatatatttctaagttttttttttttttttgtttttgtgctctgcctaagaagtgcttgagaattgtgaggagtataaaaatagtcaaagctggctgggcgcggtggctcacgcctgtaatcccagcactttgggaggctgaggcgggcggatcacgaggttaggagatggagaccatcctggctaacacagtgaaaccctgtctgtactaaagatacaaaaaactagccgggcgtggtggtaggtgcctgtagtcccagctacttgggaactcgggaggctgaggcaggagaatgacctgaacccaggaggcggagcttgcagtgagcagagattgcgccattgcactccagcctgggcgacagagtaagactgtctccaaaaaaaaaaaaaataataatcaaagctcttggatttatagtttggtccacagccttgttttgatctttcctttatcctgttttattgccatttaccacgtactgtagaaacatccctttcaactgctgataacttggaaacaagcctacaaaaataagtaatttctaactactcctaatactacctataactacccctaagcccttaccactctaacgtgacattattaaattttttattttattaacactaatattttaactacaattacagcatatgggcaatacagaatttacctaaaaggatactaatttggaacaaaaaaaatcacctttcgcacatgtatcatgtcacaaccagtttgccattgaaacaaatagaggttgcaaatattgtcagattgtcaggctgtaagaaaggatgaaattcatttcccattgcatcatcttgtggcccatggatttcaagtgccttagccaaaatcatatagctagttagcagtagagccgagactcagaaaaaaacaaagtaaaacaggcagactgaaacaaaaagtcttctaattcccagtccacatgtaaaatttgcttcatataaacaaacctaattgtaaatggcactgtagcaacaggcttctttttaacacttggattggtaaaggtcttgtttgcaacatattagaagtattatttttctctttcccccccaccccacccccaacagagtctggctctgccgcccacgccggagtgcagtggtgcagtcttggctcactgcagcctccacttcccaggttcaagcaattctcgtgcctcagcctcctgagtagttgggattacaggtgctcaccactatacccggctaatttttgtatttttagtagagatggggtttcgtcatgttggccaggttggtcttaaactcctgacctcaagtgatccacccaccttggccttccaaaatgctgggattacaggcttgagccaccaggcctttctttgttcttaggagtatagtcagactaacttctagtagttatatttctaataattgaggatgtaagtaaggatcaaatcttaaatcagtataatgcattgtcattccagagataaatcctagacccttcttggcctccttctgacataattctaatcctacagtctcagagatgctgttgtatcctgccccccaaccccatgatagtgatagtggtttttgccttgaaggaattgctttgtatttagcttttccccctctagatttctagttccttttcagtattggattggatttgagatttgattaacctagtactcaggttcagatgctcgcctctttgcaattttaacactcattcgacaataaagtcagtaaaaaacacaaaaaaaaaaaaaaaa")
                .build();
    }

    public static TranscriptModelBuilder builder() {
        return new TranscriptModelBuilder();
    }

    /**
     * Class for semi-manually building transcript models using the UCSC Table Browser.
     */
    public static class TranscriptModelBuilder {

        private GeneIdentifier geneIdentifier = null;
        private String knownGeneLine = null;
        private String mRnaSequence = null;

        public TranscriptModelBuilder geneIdentifier(GeneIdentifier geneIdentifier) {
            this.geneIdentifier = geneIdentifier;
            return this;
        }

        /**
         * http://genome-euro.ucsc.edu/cgi-bin/hgTables select output from
         * Group: Genes and Gene Predictions
         * Table: knownGene
         * paste an identifier into the 'identifiers' then press 'get output'
         * copy the single line for the accession in here.
         * @param knownGeneLine
         * @return
         */
        public TranscriptModelBuilder knownGeneLine(String knownGeneLine) {
            this.knownGeneLine = knownGeneLine;
            return this;
        }

        /**
         * http://genome-euro.ucsc.edu/cgi-bin/hgTables select output from
         * Group: Genes and Gene Predictions
         * Table: knownGeneTxMrna
         * paste an identifier into the 'identifiers' then press 'get output'
         * copy the single line for the accession in here.
         * @param mRnaSequence
         * @return
         */
        public TranscriptModelBuilder mRnaSequence(String mRnaSequence) {
            this.mRnaSequence = mRnaSequence;
            return this;
        }

        public TranscriptModel build() {
            Objects.requireNonNull(geneIdentifier);
            Objects.requireNonNull(knownGeneLine);
            Objects.requireNonNull(mRnaSequence);
            return buildTranscriptModelFromUcscKnownGene(geneIdentifier, knownGeneLine, mRnaSequence);
        }

        private TranscriptModel buildTranscriptModelFromUcscKnownGene(GeneIdentifier geneIdentifier, String knownGeneLine, String mRnaSequence) {
            String[] fields = knownGeneLine.split("\t");
            // these lines are from the UCSC knownGenes.txt.gz
            String accession = fields[0];
            int chr = parseChromosome(fields[1]);
            Strand strand = parseStrand(fields[2].charAt(0));
            int txStart = Integer.parseInt(fields[3]);
            int txEnd = Integer.parseInt(fields[4]);
            int cdsStart = Integer.parseInt(fields[5]);
            int cdsEnd = Integer.parseInt(fields[6]);
            int[][] exons = parseExons(fields[7], fields[8], fields[9]);

            Map<String, String> altGeneIds = buildAltGeneIdentifiers(geneIdentifier);

            GeneTranscriptModelBuilder geneTranscriptModelBuilder = new GeneTranscriptModelBuilder(geneIdentifier.getGeneSymbol(), geneIdentifier.getEntrezId(), accession, chr, strand, mRnaSequence)
                    .buildTxRegion(txStart, txEnd)
                    .buildCdsRegion(cdsStart, cdsEnd)
                    .addAltGeneIds(altGeneIds);

            for (int[] exon : exons) {
                geneTranscriptModelBuilder.addExon(exon[0], exon[1]);
            }

            return geneTranscriptModelBuilder.build();
        }

        private Integer parseChromosome(String field) {
            return REF_DICT.getContigNameToID().get(field.substring(3));
        }

        private Strand parseStrand(char strand) {
            return (strand == '-') ? REV : FWD;
        }

        private int[][] parseExons(String count, String startPositions, String endPositions ) {
            int exonCount = Integer.parseInt(count);
            String[] startFields = startPositions.split(",");
            String[] endFields = endPositions.split(",");

            int[][] exons = new int[exonCount][2];

            for (int i = 0; i < exonCount; ++i) {
                exons[i][0] = Integer.parseInt(startFields[i]);
                exons[i][1] = Integer.parseInt(endFields[i]);
            }

            return exons;
        }

        private Map<String,String> buildAltGeneIdentifiers(GeneIdentifier geneIdentifier) {
            Map<String, String> altGeneIds = new HashMap<>();
            altGeneIds.put(AltGeneIDType.HGNC_ID.toString(), geneIdentifier.getHgncId());
            altGeneIds.put(AltGeneIDType.HGNC_SYMBOL.toString(), geneIdentifier.getHgncSymbol());
            altGeneIds.put(AltGeneIDType.ENTREZ_ID.toString(), geneIdentifier.getEntrezId());
            altGeneIds.put(AltGeneIDType.ENSEMBL_GENE_ID.toString(), geneIdentifier.getEnsemblId());
            altGeneIds.put(AltGeneIDType.UCSC_ID.toString(), geneIdentifier.getUcscId());
            return altGeneIds;
        }
    }
}
