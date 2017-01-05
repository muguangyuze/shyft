from ._hbv_stack import *
# Fix up types that we need attached to the model
HbvStateVector.push_back = lambda self, x: self.append(x)
HbvStateVector.size = lambda self: len(self)

HbvModel.cell_t = HbvCellAll
HbvParameter.map_t = HbvParameterMap
HbvModel.parameter_t = HbvParameter
HbvModel.state_t = HbvState
HbvModel.state_with_id_t = HbvStateWithId
HbvModel.state = property(lambda self:HbvCellAllStateHandler(self.get_cells()))
HbvModel.statistics = property(lambda self: HbvCellAllStatistics(self.get_cells()))
HbvModel.hbv_snow_state = property(lambda self: HbvCellHBVSnowStateStatistics(self.get_cells()))
HbvModel.hbv_snow_response = property(lambda self: HbvCellHBVSnowResponseStatistics(self.get_cells()))
HbvModel.priestley_taylor_response = property(lambda self: HbvCellPriestleyTaylorResponseStatistics(self.get_cells()))
HbvModel.hbv_actual_evaptranspiration_response=property(lambda self: HbvCellHbvActualEvapotranspirationResponseStatistics(self.get_cells()))
HbvModel.soil_state = property(lambda self: HbvCellSoilStateStatistics(self.get_cells()))
HbvModel.tank_state = property(lambda self: HbvCellTankStateStatistics(self.get_cells()))

HbvOptModel.cell_t = HbvCellOpt
HbvOptModel.parameter_t = HbvParameter
HbvOptModel.state_t = HbvState
HbvOptModel.state_with_id_t=HbvStateWithId
HbvOptModel.state = property(lambda self:HbvCellOptStateHandler(self.get_cells()))
HbvOptModel.statistics = property(lambda self:HbvCellOptStatistics(self.get_cells()))
HbvOptModel.optimizer_t = HbvOptimizer

HbvCellAll.vector_t = HbvCellAllVector
HbvCellOpt.vector_t = HbvCellOptVector
HbvState.vector_t = HbvStateVector
HbvState.serializer_t= HbvStateIo

#decorate StateWithId for serialization support
def serialize_to_bytes(state_with_id_vector):
    if not isinstance(state_with_id_vector,HbvStateWithIdVector):
        raise RuntimeError("supplied argument must be of type HbvStateWithIdVector")
    return serialize(state_with_id_vector)

HbvStateWithIdVector.serialize_to_bytes = lambda self: serialize_to_bytes(self)

def deserialize_from_bytes(bytes):
    if not isinstance(bytes,ByteVector):
        raise RuntimeError("Supplied type must be a ByteVector, as created from serialize_to_bytes")
    states=HbvStateWithIdVector()
    deserialize(bytes,states)
    return states