//package org.monarchinitiative.exomiser.core.genome.dao.serialisers;
//
//import org.h2.mvstore.MVMap;
//import org.h2.mvstore.MVStore;
//import org.monarchinitiative.exomiser.core.proto.AlleleProto;
//
//public class MVFactory {
//
//    // Predefined mock data
//    private static final AlleleProto.AlleleKey positionStartMinus3 = AlleleProto.AlleleKey.newBuilder()
//            .setChr(1)
//            .setPosition(1227)
//            .setAlt("A")
//            .setRef("G")
//            .build();
//
//    private static final AlleleProto.AlleleKey positionEndPlus3 = AlleleProto.AlleleKey.newBuilder()
//            .setChr(1)
//            .setPosition(1233)
//            .setAlt("C")
//            .setRef("T")
//            .build();
//
//    public static MVStore createPopulatedMVStore() {
//        MVStore mvStore = new MVStore.Builder().open();
//        populateVariantMap(mvStore);
//        return mvStore;
//    }
//
//    private static void populateVariantMap(MVStore mvStore) {
//        MVMap<AlleleProto.AlleleKey, AlleleProto.ClinVar> variantMap = mvStore.openMap("clinvar");
//
//        variantMap.put(serializeKey(positionStartMinus3), mockClinVarAnnotationFor(positionStartMinus3));
//        variantMap.put(serializeKey(positionEndPlus3), mockClinVarAnnotationFor(positionEndPlus3));
//    }
//
//    private static AlleleProto.AlleleKey serializeKey(AlleleProto.AlleleKey key) {
//        return key;
//    }
//
//    private static AlleleProto.ClinVar mockClinVarAnnotationFor(AlleleProto.AlleleKey key) {
//        // Placeholder method to create a mock ClinVarAnnotation for a given AlleleKey.
//        return new AlleleProto.ClinVar("mockDataFor" + key.toString());
//    }
//}
//
