
#define interrupt_pin 2 // connect DAQ output
#define analog_pin 5 // connect DAQ output (same as interrupt_pin)
#define led_pin 13
#define test_pin 5

#define stepper_dir_pin 10 // connect to DIR terminal
#define stepper_step_pin 9 // connect to STEP terminal

#define step_delay 300

float A; // value of analog pin 0-1023
float V; // value of analog pin 0-5V
int count = 0;

void setup() {
  Serial.begin(9600);
  
  // Stepper control pins
  pinMode(stepper_dir_pin, OUTPUT);
  pinMode(stepper_step_pin, OUTPUT);
  pinMode(led_pin, OUTPUT);
  pinMode(test_pin, OUTPUT);
  digitalWrite(stepper_dir_pin, HIGH); // CW
  digitalWrite(led_pin, LOW); // off
  digitalWrite(test_pin, LOW); // off

  // DAQ read interrupt pin
  attachInterrupt(digitalPinToInterrupt(interrupt_pin), step_motor, RISING);
  Serial.print("Begin\n");
}

void loop() {
  // just wait for interrupts
  // digitalWrite(stepper_step_pin, HIGH);
  // delayMicroseconds(step_delay);
  // digitalWrite(stepper_step_pin, LOW);
  // delayMicroseconds(10);
  Serial.print(count);
  Serial.print("\n");
  delay(50);

  // digitalWrite(test_pin, HIGH);
  // delay(10);
  // digitalWrite(test_pin, LOW);
  // delay(10);
  
}

// Step the motor
void step(int delay) {
  //digitalWrite(stepper_step_pin, HIGH);
  //delayMicroseconds(delay);
  //digitalWrite(stepper_step_pin, LOW);
  //delayMicroseconds(10);
}

// Step the motor in correct direction based on value of analog pin
void step_motor() {
  //Serial.print(count);
  count++;
  //digitalWrite(led_pin, HIGH);
  //delay(1000);
  //digitalWrite(led_pin, LOW);
  
  // //step(step_delay);
  // digitalWrite(stepper_step_pin, HIGH);
  // delayMicroseconds(step_delay);
  // digitalWrite(stepper_step_pin, LOW);
  // delayMicroseconds(10);

  
  // A = analogRead(A5);  // read the input pin
  // V = map(A, 0, 1023, 0, 5); // map to voltage
  
  // Serial.println(count);

  // if (V > 3.5) { // step CW
  //   digitalWrite(stepper_dir_pin, HIGH); // CW
  //   step(step_delay);
  //   //Serial.print('CW: ');
  //   //Serial.println(V);
  // }
  // else if (V < 1.5) { // step CCW
  //   digitalWrite(stepper_dir_pin, LOW); // CCW
  //   step(step_delay);
  //   //Serial.print('CW: ');
  //   //Serial.println(V);
  // }    
  // else { // no step
  //  Serial.println('error');
  // }
}