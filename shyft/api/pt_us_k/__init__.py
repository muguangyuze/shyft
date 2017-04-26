from ._pt_us_k import *
from .. import ByteVector
# Fix up types that we need attached to the model
PTUSKStateVector.push_back = lambda self, x: self.append(x)
PTUSKStateVector.size = lambda self: len(self)


PTUSKModel.cell_t = PTUSKCellAll
PTUSKParameter.map_t = PTUSKParameterMap
PTUSKModel.parameter_t = PTUSKParameter
PTUSKModel.state_t = PTUSKState
PTUSKModel.state_with_id_t = PTUSKStateWithId
PTUSKModel.state = property(lambda self:PTUSKCellAllStateHandler(self.get_cells()))
PTUSKModel.statistics = property(lambda self: PTUSKCellAllStatistics(self.get_cells()))

PTUSKModel.universal_snow_state = property(lambda self: PTUSKCellGammaSnowStateStatistics(self.get_cells()))
PTUSKModel.universal_snow_response = property(lambda self: PTUSKCellGammaSnowResponseStatistics(self.get_cells()))
PTUSKModel.priestley_taylor_response = property(lambda self: PTUSKCellPriestleyTaylorResponseStatistics(self.get_cells()))
PTUSKModel.actual_evaptranspiration_response=property(lambda self: PTUSKCellActualEvapotranspirationResponseStatistics(self.get_cells()))
PTUSKModel.kirchner_state = property(lambda self: PTUSKCellKirchnerStateStatistics(self.get_cells()))

PTUSKOptModel.cell_t = PTUSKCellOpt
PTUSKOptModel.parameter_t = PTUSKParameter
PTUSKOptModel.state_t = PTUSKState
PTUSKOptModel.state_with_id_t = PTUSKStateWithId
PTUSKOptModel.state = property(lambda self:PTUSKCellOptStateHandler(self.get_cells()))
PTUSKOptModel.statistics = property(lambda self:PTUSKCellOptStatistics(self.get_cells()))

PTUSKOptModel.optimizer_t = PTUSKOptimizer
PTUSKOptModel.full_model_t =PTUSKModel
PTUSKModel.opt_model_t =PTUSKOptModel
PTUSKModel.create_opt_model_clone = lambda self: create_opt_model_clone(self)
PTUSKModel.create_opt_model_clone.__doc__ = create_opt_model_clone.__doc__
PTUSKOptModel.create_full_model_clone = lambda self: create_full_model_clone(self)
PTUSKOptModel.create_full_model_clone.__doc__ = create_full_model_clone.__doc__


PTUSKCellAll.vector_t = PTUSKCellAllVector
PTUSKCellOpt.vector_t = PTUSKCellOptVector
PTUSKState.vector_t = PTUSKStateVector
PTUSKState.serializer_t= PTUSKStateIo

#decorate StateWithId for serialization support
def serialize_to_bytes(state_with_id_vector):
    if not isinstance(state_with_id_vector,PTUSKStateWithIdVector):
        raise RuntimeError("supplied argument must be of type PTUSKStateWithIdVector")
    return serialize(state_with_id_vector)

PTUSKStateWithIdVector.serialize_to_bytes = lambda self: serialize_to_bytes(self)

def deserialize_from_bytes(bytes):
    if not isinstance(bytes,ByteVector):
        raise RuntimeError("Supplied type must be a ByteVector, as created from serialize_to_bytes")
    states=PTUSKStateWithIdVector()
    deserialize(bytes,states)
    return states