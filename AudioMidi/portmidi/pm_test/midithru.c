/* midithru.c -- example program implementing background thru processing */

/* suppose you want low-latency midi-thru processing, but your
   application wants to take advantage of the input buffer and
   timestamped data so that it does not have to operate with very low
   latency.

   This program illustrates how to use a timer callback from PortTime
   to implement a low-latency process that handles midi thru,
   including correctly merging midi data from the application with
   midi data from the input port.

   The main application, which runs in the main program thread, will
   use an interface similar to that of PortMidi, but since PortMidi
   does not allow concurrent threads to share access to a stream, the
   application will call private methods that transfer MIDI messages
   to and from the timer thread using lock-free queues. All PortMidi
   API calls are made from the timer thread.
 */

/* DESIGN

All setup will be done by the main thread. Then, all direct access to 
PortMidi will be handed off to the timer callback thread.

After this hand-off, the main thread will get/send messages via a queue.

The goal is to send incoming messages to the midi output while merging
any midi data generated by the application. Sysex is a problem here
because you cannot insert (merge) a midi message while a sysex is in
progress. There are at least three ways to implement midi thru with 
sysex messages:

1) Turn them off. If your application does not need them, turn them off
   with Pm_SetFilter(midi_in, PM_FILT_ACTIVE | PM_FILT_SYSEX). You will
   not receive sysex (or active sensing messages), so you will not have 
   to handle them.

2) Make them atomic. As you receive sysex messages, copy the data into
   a (big) buffer. Ideally, expand the buffer as needed -- sysex messages
   do not have any maximum length. Even more ideally, use a list structure
   and real-time memory allocation to avoid latency in the timer thread.
   When a full sysex message is received, send it to the midi output all
   at once.

3) Process sysex incrementally. Send sysex data to midi output as it
   arrives. Block any non-real-time messages from the application until
   the sysex message completes. There is the risk that an incomplete
   sysex message will block messages forever, so implement a 5-second
   timeout: if no sysex data is seen for 5 seconds, release the block,
   possibly losing the rest of the sysex message. 

   Application messages must be processed similarly: once started, a
   sysex message will block MIDI THRU processing. We will assume that
   the application will not abort a sysex message, so timeouts are not
   necessary here.

This code implements (3).

Latency is also an issue. PortMidi requires timestamps to be in 
non-decreasing order. Since we'll be operating with a low-latency
timer thread, we can just set the latency to zero meaning timestamps
are ignored by PortMidi. This will allow thru to go through with
minimal latency. The application, however, needs to use timestamps
because we assume it is high latency (the whole purpose of this
example is to illustrate how to get low-latency thru with a high-latency
application.) So the callback thread will implement midi timing by
observing timestamps. The current timestamp will be available in the
global variable current_timestamp.

*/


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "portmidi.h"
#include "pmutil.h"
#include "porttime.h"

#define MIDI_SYSEX 0xf0
#define MIDI_EOX 0xf7
#define STRING_MAX 80 /* used for console input */

/* active is set true when midi processing should start, must be
 * volatile to force thread to check for updates by other thread */
int active = FALSE;
/* process_midi_exit_flag is set when the timer thread shuts down;
 * must be volatile so it is re-read in the while loop that waits on it */
volatile int process_midi_exit_flag;

PmStream *midi_in;
PmStream *midi_out;

/* shared queues */
#define IN_QUEUE_SIZE 1024
#define OUT_QUEUE_SIZE 1024
PmQueue *in_queue;
PmQueue *out_queue;
/* this is volatile because it is set in the process_midi callback and
 * the main thread reads it to sense elapsed time. Without volatile, the
 * optimizer can put it in a register and not see the updates.
 */
volatile PmTimestamp current_timestamp = 0;
int thru_sysex_in_progress = FALSE;
int app_sysex_in_progress = FALSE;
PmTimestamp last_timestamp = 0;


static void prompt_and_exit(void)
{
    printf("type ENTER...");
    while (getchar() != '\n') ;
    /* this will clean up open ports: */
    exit(-1);
}


static PmError checkerror(PmError err)
{
    if (err == pmHostError) {
        /* it seems pointless to allocate memory and copy the string,
         * so I will do the work of Pm_GetHostErrorText directly
         */
        char errmsg[80];
        Pm_GetHostErrorText(errmsg, 80);
        printf("PortMidi found host error...\n  %s\n", errmsg);
        prompt_and_exit();
    } else if (err < 0) {
        printf("PortMidi call failed...\n  %s\n", Pm_GetErrorText(err));
        prompt_and_exit();
    }
    return err;
}


/* time proc parameter for Pm_MidiOpen */
PmTimestamp midithru_time_proc(void *info)
{
    return current_timestamp;
}


/* timer interrupt for processing midi data.
   Incoming data is delivered to main program via in_queue.
   Outgoing data from main program is delivered via out_queue.
   Incoming data from midi_in is copied with low latency to  midi_out.
   Sysex messages from either source block messages from the other.
 */
void process_midi(PtTimestamp timestamp, void *userData)
{
    PmError result;
    PmEvent buffer; /* just one message at a time */

    current_timestamp++; /* update every millisecond */

    /* do nothing until initialization completes */
    if (!active) {
        /* this flag signals that no more midi processing will be done */
        process_midi_exit_flag = TRUE;
        return;
    }

    /* see if there is any midi input to process */
    if (!app_sysex_in_progress) {
        do {
            result = Pm_Poll(midi_in);
            if (result) {
                int status;
                PmError rslt = Pm_Read(midi_in, &buffer, 1);
                if (rslt == pmBufferOverflow) 
                    continue;
                assert(rslt == 1);

                /* record timestamp of most recent data */
                last_timestamp = current_timestamp;

                /* the data might be the end of a sysex message that
                   has timed out, in which case we must ignore it.
                   It's a continuation of a sysex message if status
                   is actually a data byte (high-order bit is zero). */
                status = Pm_MessageStatus(buffer.message);
                if (((status & 0x80) == 0) && !thru_sysex_in_progress) {
                    continue; /* ignore this data */
                }

                /* implement midi thru */
                /* note that you could output to multiple ports or do other
                   processing here if you wanted
                 */
                /* printf("thru: %x\n", buffer.message); */
                Pm_Write(midi_out, &buffer, 1);

                /* send the message to the application */
                /* you might want to filter clock or active sense messages here
                   to avoid sending a bunch of junk to the application even if
                   you want to send it to MIDI THRU
                 */
                Pm_Enqueue(in_queue, &buffer);

                /* sysex processing */
                if (status == MIDI_SYSEX) thru_sysex_in_progress = TRUE;
                else if ((status & 0xF8) != 0xF8) {
                    /* not MIDI_SYSEX and not real-time, so */
                    thru_sysex_in_progress = FALSE;
                }
                if (thru_sysex_in_progress && /* look for EOX */
                    (((buffer.message & 0xFF) == MIDI_EOX) ||
                     (((buffer.message >> 8) & 0xFF) == MIDI_EOX) ||
                     (((buffer.message >> 16) & 0xFF) == MIDI_EOX) ||
                     (((buffer.message >> 24) & 0xFF) == MIDI_EOX))) {
                    thru_sysex_in_progress = FALSE;
                }
            }
        } while (result);
    }


    /* see if there is application midi data to process */
    while (!Pm_QueueEmpty(out_queue)) {
        /* see if it is time to output the next message */
        PmEvent *next = (PmEvent *) Pm_QueuePeek(out_queue);
        assert(next); /* must be non-null because queue is not empty */
        if (next->timestamp <= current_timestamp) {
            /* time to send a message, first make sure it's not blocked */
            int status = Pm_MessageStatus(next->message);
            if ((status & 0xF8) == 0xF8) {
                ; /* real-time messages are not blocked */
            } else if (thru_sysex_in_progress) {
                /* maybe sysex has timed out (output becomes unblocked) */
                if (last_timestamp + 5000 < current_timestamp) {
                    thru_sysex_in_progress = FALSE;
                } else break; /* output is blocked, so exit loop */
            }
            Pm_Dequeue(out_queue, &buffer);
            Pm_Write(midi_out, &buffer, 1);

            /* inspect message to update app_sysex_in_progress */
            if (status == MIDI_SYSEX) app_sysex_in_progress = TRUE;
            else if ((status & 0xF8) != 0xF8) {
                /* not MIDI_SYSEX and not real-time, so */
                app_sysex_in_progress = FALSE;
            }
            if (app_sysex_in_progress && /* look for EOX */
                (((buffer.message & 0xFF) == MIDI_EOX) ||
                 (((buffer.message >> 8) & 0xFF) == MIDI_EOX) ||
                 (((buffer.message >> 16) & 0xFF) == MIDI_EOX) ||
                 (((buffer.message >> 24) & 0xFF) == MIDI_EOX))) {
                app_sysex_in_progress = FALSE;
            }
        } else break; /* wait until indicated timestamp */
    }
}


void exit_with_message(char *msg)
{
#define STRING_MAX 80
    printf("%s\nType ENTER...", msg);
    while (getchar() != '\n') ;
    exit(1);
}


void initialize(int input, int output, int virtual)
/* set up midi processing thread and open midi streams */
{
    /* note that it is safe to call PortMidi from the main thread for
       initialization and opening devices. You should not make any
       calls to PortMidi from this thread once the midi thread begins.
       to make PortMidi calls.
     */

    /* note that this routine provides minimal error checking. If
       you use the PortMidi library compiled with PM_CHECK_ERRORS,
       then error messages will be printed and the program will exit
       if an error is encountered. Otherwise, you should add some
       error checking to this code.
     */

    const PmDeviceInfo *info;

    /* make the message queues */
    in_queue = Pm_QueueCreate(IN_QUEUE_SIZE, sizeof(PmEvent));
    assert(in_queue != NULL);
    out_queue = Pm_QueueCreate(OUT_QUEUE_SIZE, sizeof(PmEvent));
    assert(out_queue != NULL);

    /* always start the timer before you start midi */
    Pt_Start(1, &process_midi, 0); /* start a timer with millisecond accuracy */
    /* the timer will call our function, process_midi() every millisecond */
    
    Pm_Initialize();

    if (output < 0) {
        if (!virtual) {
            output = Pm_GetDefaultOutputDeviceID();
        }
    }
    if (output >= 0) {
        info = Pm_GetDeviceInfo(output);
        if (info == NULL) {
            printf("Could not open default output device (%d).", output);
            exit_with_message("");
        }

        printf("Opening output device %s %s\n", info->interf, info->name);

        /* use zero latency because we want output to be immediate */
        Pm_OpenOutput(&midi_out,
                      output,
                      NULL /* driver info */,
                      OUT_QUEUE_SIZE,
                      &midithru_time_proc,
                      NULL /* time info */,
                      0 /* Latency */);
    } else { /* send to virtual port */
        int id;
        printf("Opening virtual output device \"midithru\"\n");
        id = Pm_CreateVirtualOutput("midithru", NULL, NULL);
        if (id < 0) checkerror(id);  /* error reporting */
        checkerror(Pm_OpenOutput(&midi_out, id, NULL, OUT_QUEUE_SIZE,
                                 &midithru_time_proc, NULL, 0));
    }
    if (input < 0) {
        if (!virtual) {
            input = Pm_GetDefaultInputDeviceID();
        }
    }
    if (input >= 0) {
        info = Pm_GetDeviceInfo(input);
        if (info == NULL) {
            printf("Could not open default input device (%d).", input);
            exit_with_message("");
        }
        
        printf("Opening input device %s %s\n", info->interf, info->name);
        Pm_OpenInput(&midi_in,
                     input,
                     NULL /* driver info */,
                     0 /* use default input size */,
                     &midithru_time_proc,
                     NULL /* time info */);
    } else { /* receive from virtual port */
        int id;
        printf("Opening virtual input device \"midithru\"\n");
        id = Pm_CreateVirtualInput("midithru", NULL, NULL);
        if (id < 0) checkerror(id);  /* error reporting */
        checkerror(Pm_OpenInput(&midi_in, id, NULL, 0,
                                &midithru_time_proc, NULL));
    }
    /* Note: if you set a filter here, then this will filter what goes
       to the MIDI THRU port. You may not want to do this.
     */
    Pm_SetFilter(midi_in, PM_FILT_ACTIVE | PM_FILT_CLOCK);

    active = TRUE; /* enable processing in the midi thread -- yes, this
                      is a shared variable without synchronization, but
                      this simple assignment is safe */

}


void finalize()
{
    /* the timer thread could be in the middle of accessing PortMidi stuff */
    /* to detect that it is done, we first clear process_midi_exit_flag and
       then wait for the timer thread to set it
     */
    process_midi_exit_flag = FALSE;
    active = FALSE;
    /* busy wait for flag from timer thread that it is done */
    while (!process_midi_exit_flag) ;
    /* at this point, midi thread is inactive and we need to shut down
     * the midi input and output
     */
    Pt_Stop(); /* stop the timer */
    Pm_QueueDestroy(in_queue);
    Pm_QueueDestroy(out_queue);

    Pm_Close(midi_in);
    Pm_Close(midi_out);

    Pm_Terminate();    
}


int main(int argc, char *argv[])
{
    PmTimestamp last_time = 0;
    PmEvent buffer;
    int i;
    int input = -1, output = -1;
    int virtual = FALSE;
    int delay_enable = TRUE;

    printf("Usage: midithru [-i input] [-o output] [-v] [-n]\n"
           "where input and output are portmidi device numbers\n"
           "if -v and input and/or output are not specified,\n"
           "then virtual ports are created and used instead.\n"
           "-n turns off the default MIDI delay effect.\n");
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            i++;
            input = atoi(argv[i]);
            printf("Input device number: %d\n", input);
        } else if (strcmp(argv[i], "-o") == 0) {
            i++;
            output = atoi(argv[i]);
            printf("Output device number: %d\n", output);
        } else if (strcmp(argv[i], "-v") == 0) {
            virtual = TRUE;
        } else if (strcmp(argv[i], "-n") == 0) {
            delay_enable = FALSE;
            printf("delay_effect is disabled\n");
        } else {
            return -1;
        }
    }
    printf("begin PortMidi midithru program...\n");

    initialize(input, output, virtual); /* set up and start midi processing */
	
    printf("This program will run for 60 seconds, "
           "or until you play B below middle C,\n"
           "All input is sent immediately, implementing software MIDI THRU.\n"
           "Also, all input is echoed with a 2 second delay.\n");

    while (current_timestamp < 60000) {
        /* just to make the point that this is not a low-latency process,
           spin until half a second has elapsed */
        last_time = last_time + 500;
        while (last_time > current_timestamp) ;

        /* now read data and send it after changing timestamps */
        while (Pm_Dequeue(in_queue, &buffer) == 1) {
            /* printf("timestamp %d\n", buffer.timestamp); */
            /* printf("message %x\n", buffer.message); */
            if (delay_enable) {
                buffer.timestamp = buffer.timestamp + 2000; /* delay */
                Pm_Enqueue(out_queue, &buffer);
            }
            /* play B3 to break out of loop */
            if (Pm_MessageStatus(buffer.message) == 0x90 &&
                Pm_MessageData1(buffer.message) == 59) {
                goto quit_now;
            }
        }
    }
quit_now:
    finalize();
    exit_with_message("finished PortMidi midithru program.");
    return 0; /* never executed, but keeps the compiler happy */
}