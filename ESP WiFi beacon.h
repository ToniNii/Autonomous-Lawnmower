/* Code for ESP32 to operate as wifi beacon. 
	Simply creates AP with specified SSID. 
	!UNTESTED!*/
	
#include <WiFi.h>

const char *ssid = "SSID"					//insert intended AP SSID here
const char *password = "PASSWORD"			//insert password for AP here

void setup() {

Serial.begin(115200);						//Serial monitor initiated for debugging purposes

WiFi.softAP(ssid, password);				//Starts the AP

Serial.println("AP Active, IP addr: ");		
Serail.println(WiFi.softAPIP());			//AP IP address printed to serial monitor
}