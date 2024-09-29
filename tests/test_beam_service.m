% tests/test_beam_service.m
function test_beam_service
    % Create a Gaussian beam
    beam = ParaxialBeam(633e-9, 1e-3, 0);
    
    % Pass it to the BeamService
    beamService = BeamService(beam);
    
    % Test the beam waist calculation
    z = 1;
    expected_waist = beam.beam_waist(z);
    assert(abs(beamService.calculateBeamWaist(z) - expected_waist) < 1e-9, ...
           'Beam waist calculation is incorrect.');
       
    disp('All tests passed!');
end
