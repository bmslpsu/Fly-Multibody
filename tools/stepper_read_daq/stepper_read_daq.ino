
#define daq_step_pin 2 // step pin from DAQ
// #define daq_dir_pin 3 // dir pin from DAQ

// #define stepper_dir_pin 10 // connect to DIR terminal
#define stepper_step_pin 9 // connect to STEP terminal

#define step_delay 150

#define led_pin 13
#define test_pin 5

volatile unsigned int count = 0;

void setup() {
  Serial.begin(9600);
  
  // DAQ read interrupt pin
  attachInterrupt(digitalPinToInterrupt(daq_step_pin), step, RISING );
  
  // Stepper control pins
  pinMode(stepper_step_pin, OUTPUT);

  Serial.print("Begin\n");
}

void loop() {
  // just wait for interrupts
  Serial.print(count);
  Serial.print("\n");
  delay(2000);
}

// Step the motor
void step() {
  digitalWrite(stepper_step_pin, HIGH);
  delayMicroseconds(step_delay);
  digitalWrite(stepper_step_pin, LOW);
  delayMicroseconds(10);
  count++;
}
