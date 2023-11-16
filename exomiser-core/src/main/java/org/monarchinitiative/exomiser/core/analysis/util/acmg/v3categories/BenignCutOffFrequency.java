package org.monarchinitiative.exomiser.core.analysis.util.acmg.v3categories;

import java.util.Collections;
import java.util.List;

public class BenignCutOffFrequency {

    private List<Double> benignFrequencies;
    private List<Double> pathogenicFrequencies;

    public BenignCutOffFrequency(List<Double> benignFrequencies, List<Double> pathogenicFrequencies) {
        this.benignFrequencies = benignFrequencies;
        this.pathogenicFrequencies = pathogenicFrequencies;
    }


//    find the highest frequency in the list of pathogenic variant frequencies.
//Find the Lowest Benign Frequency That Exceeds the Highest Pathogenic Frequency
// and return the first benign frequency that is greater than the highest pathogenic frequency.

    public double calculateCutoffFrequency() {
        double maxPathogenicFrequency = Collections.max(pathogenicFrequencies);

        for (Double benignFreq : benignFrequencies) {
            if (benignFreq > maxPathogenicFrequency) {
                return benignFreq;
            }
        }
        return -1;
    }
}

