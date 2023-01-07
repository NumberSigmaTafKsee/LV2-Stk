//  mutex, thread not implemented now.
%module stk 
%{
#define __OS_LINUX__
#include <Stk.h>
#include <Generator.h>
#include <Instrmnt.h>
#include <FM.h>
#include <Filter.h>
#include <Effect.h>
#include <Function.h>
#include <WvIn.h>
#include <WvOut.h>
#include <Sampler.h>
#include <ADSR.h>
#include <Asymp.h>
#include <BandedWG.h>
#include <BeeThree.h>
#include <BiQuad.h>
#include <Blit.h>
#include <BlitSaw.h>
#include <BlitSquare.h>
#include <BlowBotl.h>
#include <BlowHole.h>
#include <BowTable.h>
#include <Bowed.h>
#include <Brass.h>
#include <Chorus.h>
#include <Clarinet.h>
#include <Cubic.h>
#include <Delay.h>
#include <DelayA.h>
#include <DelayL.h>
#include <Drummer.h>
#include <Echo.h>
#include <Effect.h>
#include <Envelope.h>
#include <FMVoices.h>
#include <FileRead.h>
#include <FileWrite.h>
#include <FileWvIn.h>
#include <FileWvOut.h>
#include <FileLoop.h>
#include <Filter.h>
#include <Fir.h>
#include <Flute.h>
#include <FormSwep.h>
#include <FreeVerb.h>
#include <Granulate.h>
#include <Guitar.h>
#include <HevyMetl.h>
#include <Iir.h>
#include <JCRev.h>
#include <JetTable.h>
#include <LentPitShift.h>
#include <Mandolin.h>
#include <Mesh2D.h>
#include <Messager.h>
#include <MidiFileIn.h>
#include <Modal.h>
#include <ModalBar.h>
#include <Modulate.h>
#include <Moog.h>
#include <NRev.h>
#include <Noise.h>
#include <OnePole.h>
#include <OneZero.h>
#include <PRCRev.h>
#include <PercFlut.h>
#include <Phonemes.h>
#include <PitShift.h>
#include <Plucked.h>
#include <PoleZero.h>
#include <Recorder.h>
#include <ReedTable.h>
#include <Resonate.h>
#include <Rhodey.h>
#include <Saxofony.h>
#include <Shakers.h>
#include <Simple.h>
#include <SineWave.h>
#include <SingWave.h>
#include <Sitar.h>
#include <Sphere.h>
#include <StifKarp.h>
#include <TapDelay.h>
#include <TubeBell.h>
#include <Twang.h>
#include <TwoPole.h>
#include <TwoZero.h>
#include <Vector3D.h>
#include <VoicForm.h>
#include <Voicer.h>
#include <Whistle.h>
#include <Wurley.h>
#include <Socket.h>
#include <Mutex.h>
#include <Thread.h>
#include <InetWvIn.h>
#include <InetWvOut.h>
#include <TcpClient.h>
#include <TcpServer.h>
#include <UdpSocket.h>
#include <RtWvIn.h>
#include <RtWvOut.h>

using namespace stk;
%}

%include "std_string.i"
%include "std_except.i"

%include <Stk.h>
%import <Generator.h>
%import <Instrmnt.h>
%import <FM.h>
%import <Filter.h>
%import <Effect.h>
%import <Function.h>
%import <WvIn.h>
%import <WvOut.h>
%import <Sampler.h>

%include <ADSR.h>
%include <Asymp.h>
%include <BandedWG.h>
%include <BeeThree.h>
%include <BiQuad.h>
%include <Blit.h>
%include <BlitSaw.h>
%include <BlitSquare.h>
%include <BlowBotl.h>
%include <BlowHole.h>
%include <BowTable.h>
%include <Bowed.h>
%include <Brass.h>
%include <Chorus.h>
%include <Clarinet.h>
%include <Cubic.h>
%include <Delay.h>
%include <DelayA.h>
%include <DelayL.h>
%include <Drummer.h>
%include <Echo.h>
%include <Effect.h>
%include <Envelope.h>
%include <FMVoices.h>
%include <FileRead.h>
%include <FileWrite.h>
%include <FileWvIn.h>
%include <FileWvOut.h>
%include <FileLoop.h>
%include <Filter.h>
%include <Fir.h>
%include <Flute.h>
%include <FormSwep.h>
%include <FreeVerb.h>
%include <Granulate.h>
%include <Guitar.h>
%include <HevyMetl.h>
%include <Iir.h>
%include <JCRev.h>
%include <JetTable.h>
%include <LentPitShift.h>
%include <Mandolin.h>
%include <Mesh2D.h>
%include <Messager.h>
%include <MidiFileIn.h>
%include <Modal.h>
%include <ModalBar.h>
%include <Modulate.h>
%include <Moog.h>
%include <NRev.h>
%include <Noise.h>
%include <OnePole.h>
%include <OneZero.h>
%include <PRCRev.h>
%include <PercFlut.h>
%include <Phonemes.h>
%include <PitShift.h>
%include <Plucked.h>
%include <PoleZero.h>
%include <Recorder.h>
%include <ReedTable.h>
%include <Resonate.h>
%include <Rhodey.h>
//%include <SKINImsg.h>
//%include <SKINItbl.h>
%include <Saxofony.h>
%include <Shakers.h>
%include <Simple.h>
%include <SineWave.h>
%include <SingWave.h>
%include <Sitar.h>
//%include <Skini.h>
%include <Sphere.h>
%include <StifKarp.h>
%include <TapDelay.h>
%include <TubeBell.h>
%ignore stk::Twang::setLoopFilter;
%include <Twang.h>
%include <TwoPole.h>
%include <TwoZero.h>
%include <Vector3D.h>
%include <VoicForm.h>
%include <Voicer.h>
%include <Whistle.h>
%include <Wurley.h>
//%include <Mutex.h>
//%include <Thread.h>
%include <Socket.h>
%include <InetWvIn.h>
%include <InetWvOut.h>
%include <TcpClient.h>
%include <TcpServer.h>
%include <UdpSocket.h>
%include <RtWvIn.h>
%include <RtWvOut.h>
