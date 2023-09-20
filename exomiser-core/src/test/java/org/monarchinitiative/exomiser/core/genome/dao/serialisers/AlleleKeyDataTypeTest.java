package org.monarchinitiative.exomiser.core.genome.dao.serialisers;

import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import static org.junit.jupiter.api.Assertions.*;

class AlleleKeyDataTypeTest {

    private final  AlleleProto.AlleleKey positionStartMinus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1229)
            .setAlt("A")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionEndPlus1 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1231)
            .setAlt("C")
            .setRef("G")
            .build();

    // Variants directly at the boundaries:
    private final  AlleleProto.AlleleKey positionStartMinus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1228)
            .setAlt("A")
            .setRef("G")
            .build();

    private final  AlleleProto.AlleleKey positionEndPlus2 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1232)
            .setAlt("T")
            .setRef("A")
            .build();


    // Variants just outside the boundaries:
    private final  AlleleProto.AlleleKey positionStartMinus3 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1227)
            .setAlt("A")
            .setRef("G")
            .build();

    private final AlleleProto.AlleleKey positionEndPlus3 = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1233)
            .setAlt("C")
            .setRef("T")
            .build();

    // Variants are exactly matching on position:
    private final  AlleleProto.AlleleKey positionExactlyMatches = AlleleProto.AlleleKey.newBuilder()
            .setChr(1)
            .setPosition(1230)
            .setAlt("A")
            .setRef("G")
            .build();




    public int compare(AlleleProto.AlleleKey a, AlleleProto.AlleleKey b) {
        return AlleleKeyDataType.INSTANCE.compare(a, b);
    }

    @Test
    public void testCompare() {
        assertTrue(compare(positionStartMinus1, positionEndPlus1) < 0); // Expects positionStartMinus1 to be less than positionEndPlus1
        assertTrue(compare(positionStartMinus2, positionStartMinus1) < 0); // Expects positionStartMinus2 to be less than positionStartMinus1
        assertTrue(compare(positionEndPlus1, positionEndPlus2) < 0); // Expects positionEndPlus1 to be less than positionEndPlus2
        assertTrue(compare(positionStartMinus3, positionStartMinus2) < 0); // Expects positionStartMinus3 to be less than positionStartMinus2
        assertTrue(compare(positionEndPlus2, positionEndPlus3) < 0); // Expects positionEndPlus2 to be less than positionEndPlus3
        assertTrue(compare(positionExactlyMatches, positionEndPlus1) < 0); // Expects positionExactlyMatches to be less than positionEndPlus1
        assertTrue(compare(positionExactlyMatches, positionExactlyMatches) == 0); // Expects positionExactlyMatches to be equal to itself
    }

}